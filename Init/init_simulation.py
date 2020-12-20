# Initiate the crystal based on the chosen variables
# This will eventually become "Initiate the system" => system depends on user's choice
import os
import shutil
from .init_values import *
from tkinter import *
import Calculations.calculations as calc
from Optimade.optimade import translate_to_optimade, concatenateOptimadeDataFiles
import Json.write_json as write_json
from asap3 import Trajectory
from ase.gui import *


def simulation(EMT_Check,openKIM_Check, Lennard_Jones_Check, LJ_epsilon,
                        LJ_sigma, LJ_cutoff, Verlocity_Verlet_Check, 
                        Langevin_Check, Langevin_friction, time_step, KIM_potential,
                        ASE, Symbol, Materials_project,API_Key,Criteria_list, 
                        Vacancy, Impurity, Impurity_ele_list,
                        Temperature, Steps, Interval,Size_X, Size_Y, Size_Z,
                        PBC, Bravais_lattice,Directions,Miller,
                        lc_a,lc_b,lc_c,lc_alpha,lc_beta,lc_gamma,run_Optimade,Optimade_name):
    
    """ Function that looks if the user wants to run ASE or Materials_project 
        Checks if the simulation is going to add impurites or not
        Impurites doesnt work with openKIM however
    """
    if (ASE == True) and (Materials_project == False):
        if Impurity == True:
            for Impurity_ele in Impurity_ele_list: 
                atoms = init(EMT_Check, openKIM_Check, Lennard_Jones_Check, LJ_epsilon,
                            LJ_sigma, LJ_cutoff,Verlocity_Verlet_Check, KIM_potential,Symbol,
                            Vacancy, Impurity, Impurity_ele, Temperature,
                            Size_X,Size_Y,Size_Z,PBC,Bravais_lattice, Directions,Miller,
                            lc_a,lc_b,lc_c,lc_alpha,lc_beta,lc_gamma)

        else:
            atoms = init(EMT_Check, openKIM_Check, Lennard_Jones_Check, LJ_epsilon,
                            LJ_sigma, LJ_cutoff,Verlocity_Verlet_Check, KIM_potential,Symbol,
                            Vacancy, Impurity, Impurity_ele_list, Temperature,
                            Size_X,Size_Y,Size_Z, PBC, Bravais_lattice, Directions,Miller,
                            lc_a,lc_b,lc_c,lc_alpha,lc_beta,lc_gamma)
    elif (Materials_project == True) and (ASE == False):
        if Impurity == True:
            for Impurity_ele in Impurity_ele_list:
                atoms = init_MP(EMT_Check,openKIM_Check,Lennard_Jones_Check, LJ_epsilon,
                                LJ_sigma, LJ_cutoff,Verlocity_Verlet_Check,KIM_potential,Criteria_list,
                                Vacancy, Impurity, Impurity_ele, Temperature,
                                Size_X,Size_Y,Size_Z,API_Key,PBC)

        else:
            atoms = init_MP(EMT_Check,openKIM_Check,Lennard_Jones_Check, LJ_epsilon,
                                LJ_sigma, LJ_cutoff,Verlocity_Verlet_Check,KIM_potential,Criteria_list,
                                Vacancy, Impurity, Impurity_ele_list, Temperature,
                                Size_X,Size_Y,Size_Z,API_Key,PBC)
    else:
        raise Exception("ASE=Materials_Project. Both cannot be true/false at the same time!")
    
    count = 0 # To get a unique name for all the files

    for atomobj in atoms:
        if (Verlocity_Verlet_Check == True) and (Langevin_Check == False):
            # We want to run MD with constant energy using the VelocityVerlet algorithm.
            dyn = VelocityVerlet(atomobj, time_step*units.fs)
            
        elif (Verlocity_Verlet_Check == False) and (Langevin_Check == True):
            dyn = Langevin(atomobj, time_step*units.fs, units.kB*Temperature, Langevin_friction)
        else:
            raise Exception("Velocity_Verlet=Langevin. Both cannot be true/false at the same time!")
        #Creates a unique name for every simulation run 
        trajFileName = atomobj.get_chemical_formula() +"_run" + str(count) +  '_.traj'
        traj = Trajectory(trajFileName, "w", atomobj)
        dyn.attach(traj.write, Interval)
        dyn.run(Steps)
        traj.close()
        count = count + 1  

        traj = Trajectory(trajFileName)

        latticeConstant_a = calc.calc_lattice_constant_fcc_cubic(Symbol, EMT())

        # if EMT_Check == True:
        #     latticeConstant_a = calc.calc_lattice_constant_fcc_cubic(Symbol, EMT())
        # elif openKIM_Check == True:
        #     potential = checkKIMpotential(KIM_potential)
        #     latticeConstant_a = calc.calc_lattice_constant_fcc_cubic(Symbol, KIM(potential))
        # elif Lennard_Jones_Check == True:
        #     latticeConstant_a = calc.calc_lattice_constant_fcc_cubic(Symbol, LennardJones(list(dict.fromkeys(atomobj.get_atomic_numbers())), LJ_epsilon, LJ_sigma, rCut=LJ_cutoff, modified=True))
        
        print("Lattice constant a:", latticeConstant_a) 
        
        eq_index = calc.eq_test(atomobj, traj)
        
        if eq_index != 0:#If-statement that checks if we ever reached equilibrium. Returns a message if the traj-file is empty, otherwise does calculations.
            meansSquareDisplacement = calc.MSD_calc(atomobj, traj, -1, eq_index)
            print("MSD = ", meansSquareDisplacement, "[Å²]")
            
            selfDiffusionCoffecient = calc.Self_diffuse(meansSquareDisplacement, (len(traj) - eq_index), Interval, time_step)
            print("D = ", selfDiffusionCoffecient, "[Å²/fs]")
            
            lindemann = calc.Lindemann(traj, meansSquareDisplacement)
            print("Lindeman Criterion = ", lindemann)
            
            if Verlocity_Verlet_Check == True:    
                specificHeatCapacity = calc.Heat_Capcity_NVE(atomobj, traj, eq_index)
                print("C_v = ", specificHeatCapacity, "[J/K*Kg]")
            else:
                specificHeatCapacity = calc.Heat_Capcity_NVT(atomobj, traj, eq_index)
                print("C_v = ", specificHeatCapacity, "[J/K*Kg]")
            
            internalTemperature = calc.internal_temperature(atomobj, traj, eq_index)
            print("Internal temperature:", internalTemperature, "[K]")
            
            cohesiveEnergy = calc.cohesive_energy(atomobj, traj, eq_index)
            print("Cohesive energy:", cohesiveEnergy, "[eV/atom]")
            
            internalPressure = calc.calc_internal_pressure(atomobj, traj, eq_index, Size_X * Size_Y * Size_Z)
            print("Internal Pressure:", internalPressure, "[eV / Å^3]")
            
            debyeTemperature = calc.debye_temperature(traj,meansSquareDisplacement,eq_index)
            print("Debye temperature:", debyeTemperature, "[K]")

            e0, v0, B_GPa = calc.calc_bulk_modulus(atomobj)
            print('Bulk Modulus:', B_GPa, '[GPa]', '|', 'Minimum energy E =', e0, '[eV], at volume V =', v0, '[Å^3].')

            if run_Optimade == True:
                translate_to_optimade(atomobj, meansSquareDisplacement, selfDiffusionCoffecient, lindemann , specificHeatCapacity, 
                                        internalTemperature, cohesiveEnergy, internalPressure, B_GPa)

                concatenateOptimadeDataFiles(Optimade_name) ### Move me to supercomputer script later!!! 
            
            #Writes a .csv-file with time evolution of the mean square displacement
            calc.write_time_evolution_to_csv(atomobj, "Visualization/properties.csv", traj, eq_index, Interval, time_step)

            #Parameters needed for write_simulation_to_json. These if-checks might not be necessary if I just save the boolean instead of the name.

            database = ""
            potential = ""
            integrator = ""

            if EMT_Check == False:
                if KIM_potential == " ":
                    potential = "Leonnard-Jones"
                else:
                  potential = KIM_potential
            else:
                potential = "EMT"

            if Verlocity_Verlet_Check == True:
                integrator = "Velocity-Verlet"
            else:
                integrator = "Langevin"
                
            if Materials_project == True:
                database = "Materials Project"
            else:
                database = "ASE"

            #Makes a simulation-specific json-file with all the relevant input and output.  
            write_json.write_simulation_to_json(database,
                                                potential,
                                                integrator,
                                                atomobj,
                                                Temperature,
                                                meansSquareDisplacement,
                                                selfDiffusionCoffecient,
                                                lindemann,
                                                specificHeatCapacity,
                                                internalTemperature,
                                                cohesiveEnergy,
                                                internalPressure,
                                                B_GPa,
                                                latticeConstant_a,
                                                Steps*time_step)

        else:
            print("System never reached equilibrium. No calculations are possible.")

        #Moves the trajectory file to another folder after it has been used
        shutil.move(trajFileName, "Traj/" + trajFileName)
        
    atoms.clear()
