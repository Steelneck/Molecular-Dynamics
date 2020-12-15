# Initiate the crystal based on the chosen variables
# This will eventually become "Initiate the system" => system depends on user's choice
import os
import shutil
from .init_values import *
from tkinter import *
import Calculations.calculations as calc
import Json.write_json as write_json
from asap3 import Trajectory
from ase.gui import *


def simulation(EMT_Check,openKIM_Check,KIM_potential, Verlocity_Verlet_Check, Langevin_Check,
                        ASE, Materials_project,Symbol,Criteria_list, 
                        Vacancy, Impurity, Impurity_ele_list,
                        Temperature, Steps, Interval,
                        Size_X, Size_Y, Size_Z,API_Key,PBC,Directions,Miller,
                        lc_a,lc_b,lc_c,lc_alpha,lc_beta,lc_gamma):
    
    """ Function that looks if the user wants to run ASE or Materials_project 
        Checks if the simulation is going to add impurites or not
        Impurites doesnt work with openKIM however
    """
    if (ASE == True) and (Materials_project == False):
        if Impurity == True:
            for Impurity_ele in Impurity_ele_list: 
                atoms = init(EMT_Check, openKIM_Check, Verlocity_Verlet_Check, KIM_potential,Symbol,
                            Vacancy, Impurity, Impurity_ele, Temperature,
                            Size_X,Size_Y,Size_Z,PBC,Directions,Miller,
                            lc_a,lc_b,lc_c,lc_alpha,lc_beta,lc_gamma)

        else:
            atoms = init(EMT_Check, openKIM_Check, Verlocity_Verlet_Check, KIM_potential,Symbol,
                            Vacancy, Impurity, Impurity_ele_list, Temperature,
                            Size_X,Size_Y,Size_Z,PBC,Directions,Miller,
                            lc_a,lc_b,lc_c,lc_alpha,lc_beta,lc_gamma)
    elif (Materials_project == True) and (ASE == False):
        if Impurity == True:
            for Impurity_ele in Impurity_ele_list:
                atoms = init_MP(EMT_Check,openKIM_Check,Verlocity_Verlet_Check,KIM_potential,Criteria_list,
                                Vacancy, Impurity, Impurity_ele, Temperature,
                                Size_X,Size_Y,Size_Z,API_Key,PBC)

        else:
            atoms = init_MP(EMT_Check,openKIM_Check,Verlocity_Verlet_Check,KIM_potential,Criteria_list,
                                Vacancy, Impurity, Impurity_ele_list, Temperature,
                                Size_X,Size_Y,Size_Z,API_Key,PBC)
    else:
        raise Exception("ASE=Materials_Materials. Both cannot be true/false at the same time!")
    
    for atomobj in atoms:

        if (Verlocity_Verlet_Check == True) and (Langevin_Check == False):
            # We want to run MD with constant energy using the VelocityVerlet algorithm.
            dyn = VelocityVerlet(atomobj, 5*units.fs)  # 5 fs time step.
        elif (Verlocity_Verlet_Check == False) and (Langevin_Check == True):
            dyn = Langevin(atomobj, 5*units.fs, units.kB*Temperature, 0.002)
        else:
            raise Exception("Velocity_Verlet=Langevin. Both cannot be true/false at the same time!")

        #Creates a unique name for every simulation run 
        trajFileName = atomobj.get_chemical_formula() + '.traj'

        traj = Trajectory(trajFileName, "w", atomobj)
        
        dyn.attach(traj.write, Interval)
        dyn.run(Steps)
        
        traj = Trajectory(trajFileName)

        latticeConstant_a = calc.calc_lattice_constant_fcc_cubic(Symbol, EMT())
        
        eq_index = calc.eq_traj(atomobj, traj, Size_X * Size_Y * Size_Z)
        if eq_index != 0:
        #if os.path.getsize(trajFileName_eq) != 0: #If-statement that checks if we ever reached equilibrium. Returns a message if the traj-file is empty, otherwise does calculations.
            #calc.write_atom_properties(atoms, "Visualization/properties.csv", traj, eq_index)
            
            MSD = calc.MSD_calc(atomobj, traj, eq_index)
            print("MSD = ", MSD, "[Å²]")
            
            D = calc.Self_diffuse(MSD, (len(traj) - eq_index), Interval)
            print("D = ", D, "[Å²/fs]")
            
            L = calc.Lindemann(traj, MSD)
            if L>0.1:
                print(L, "Melting according to Lindemann criterion.")
            else:
                print(L, "Not melting according to Lindemann criterion.")
            
            SHC = calc.Specific_Heat(atomobj, traj, eq_index)
            print("C_p = ", SHC, "[J/K*Kg]")
            
            internalTemperature = calc.internal_temperature(atomobj, traj, eq_index)
            print("Internal temperature:", internalTemperature, "[K]")
            
            cohesiveEnergy = calc.cohesive_energy(atomobj, traj, eq_index)
            print("Cohesive energy:", cohesiveEnergy, "[eV/atom]")
            
            internalPressure = calc.calc_internal_pressure(atomobj, traj, eq_index, Size_X * Size_Y * Size_Z)
            print("Internal Pressure:", internalPressure, "[eV / Å^3]")
            
            e0, v0, B_GPa = calc.calc_bulk_modulus(atomobj)
            print('Bulk Modulus:', B_GPa, '[GPa]', '|', 'Minimum energy E =', e0, '[eV], at volume V =', v0, '[Å^3].')

            #Moves the trajectory file to another folder after it has been used
            shutil.move(trajFileName, "Traj/" + trajFileName)

            #Makes a simulation-specific json-file with all the relevant input and output.

            database = ""
            ele_symbol = ""
            potential = ""
            integrator = ""
            
            if Materials_project == True and ASE == False:
                ele_symbol = Criteria_list[0]["elements"][0]
            else:
                ele_symbol = Symbol
                
            if EMT == True:
                potential = "EMT"
            elif KIM_potential == " ":
                potential = "Leonnard-Jones"
            else:
                potential = KIM_potential

            if Verlocity_Verlet_Check == True:
                integrator = "Velocity-Verlet"
            else:
                integrator = "Langevin"
                
            if Materials_project == True:
                database = "Materials Project"
            else:
                database = "ASE"
                
            write_json.write_simulation_to_json(ele_symbol,
                                                database,
                                                potential,
                                                integrator,
                                                atomobj,
                                                Temperature,
                                                MSD,
                                                D,
                                                L,
                                                SHC,
                                                internalTemperature,
                                                cohesiveEnergy,
                                                internalPressure,
                                                B_GPa,
                                                latticeConstant_a,
                                                Steps*5)

        else:
            print("System never reached equilibrium. No calculations are possible.")
    atoms.clear()
