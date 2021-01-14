# Initiate the crystal based on the chosen variables
# This will eventually become "Initiate the system" => system depends on user's choice
import os
import shutil
from .init_values import *
from tkinter import *
import Json.write_json as write_json

from .init_functions import create_vacancy, find_crystal_center, insert_impurity
from .init_values import checkKIMpotential
import Calculations.calculations as calc
from Optimade.optimade import translate_to_optimade, concatenateOptimadeDataFiles

from ase.gui import *
from ase import units
#from ase.calculators.kim.kim import KIM

# Algorithms and calculators for the simulation
from asap3.md.velocitydistribution import MaxwellBoltzmannDistribution, Stationary
from ase.md.velocitydistribution import ZeroRotation
from asap3.md.verlet import VelocityVerlet
from asap3.md.langevin import Langevin
from asap3 import EMT
from asap3 import Trajectory
from asap3 import LennardJones
from asap3 import OpenKIMcalculator 


def simulation(EMT_Check,openKIM_Check, Lennard_Jones_Check, LJ_epsilon,
                        LJ_sigma, LJ_cutoff, Velocity_Verlet_Check, 
                        Langevin_Check, Langevin_friction, Time_step, KIM_potential,
                        ASE, Symbol, Materials_project,API_Key,Criteria_list, 
                        Vacancy, Impurity, Impurity_ele,
                        Temperature, Steps, Interval,Size_X, Size_Y, Size_Z,
                        PBC, Bravais_lattice,Directions,Miller,
                        lc_a,lc_b,lc_c,lc_alpha,lc_beta,lc_gamma,Run_Optimade,Optimade_name,Optimized_volume):
    
    """ Function that looks if the user wants to run ASE or Materials_project 
        Checks if the simulation is going to add impurites or not
        Impurites doesnt work with openKIM however
    """
    if (ASE == True) and (Materials_project == False):

        atoms = init(Symbol, Bravais_lattice, Directions, Miller,
                        lc_a,lc_b,lc_c,lc_alpha,lc_beta,lc_gamma)


    elif (Materials_project == True) and (ASE == False):
        atoms = init_MP(Criteria_list, API_Key)

    else:
        raise Exception("ASE=Materials_Project. Both cannot be true/false at the same time!")
    
    count = 0 # To get a unique name for all the files
    print(len(atoms))
    for atomobj in atoms:
        try: 
            # Run simulation with optimized volume.
            if Optimized_volume == True:

                # Save all the angles for the unit cell
                alpha = (atomobj.get_cell_lengths_and_angles())[3]
                beta = (atomobj.get_cell_lengths_and_angles())[4]
                gamma = (atomobj.get_cell_lengths_and_angles())[5]

                #Calculates the optimal length for the lattice constant 
                if EMT_Check == True:
                    latticeConstant_c = calc.calc_lattice_constant_cubic(atomobj, EMT(), alpha, beta, gamma, Size_X, Size_Y, Size_Z, PBC)
                elif openKIM_Check == True:
                    potential = checkKIMpotential(KIM_potential)
                    #latticeConstant_c = calc.calc_lattice_constant_cubic(atomobj, KIM(potential, options={"ase_neigh": True}), 
                    #                                                                         alpha, beta, gamma, Size_X, Size_Y, Size_Z, PBC)
                    latticeConstant_c = calc.calc_lattice_constant_cubic(atomobj, OpenKIMcalculator(potential), 
                                                                                            alpha, beta, gamma, Size_X, Size_Y, Size_Z, PBC)
                elif Lennard_Jones_Check == True:
                    latticeConstant_c = calc.calc_lattice_constant_cubic(atomobj, LennardJones(list(dict.fromkeys(atomobj.get_atomic_numbers())), 
                                                                                            LJ_epsilon, LJ_sigma, rCut=LJ_cutoff, modified=True), alpha, beta, gamma, Size_X, Size_Y, Size_Z, PBC)
                # Updates the cell with the optimal lattice constant
                atomobj.set_cell([latticeConstant_c, latticeConstant_c, latticeConstant_c, alpha, beta, gamma])
                print("lattice constant:", latticeConstant_c, "\n")
            else:
                # The lattice constant is set to zero if Optimized_volume == false. This is implemented so that the visualization does not crash!
                latticeConstant_c = 0 

            #Creates a supercell
            atomobj = atomobj*(Size_X,Size_Y,Size_Z)

            #Set the periodic boundary conditions
            atomobj.set_pbc(PBC)

            #places impurity in the crystal 
            if Impurity == True:
                atom_pos = find_crystal_center(atomobj) # Returns a center position in the crystal
                insert_impurity(atomobj, Impurity_ele, atom_pos) # Insert "foregin" atom in the crystal

            #Places vacancy in the crytal
            if Vacancy == True:
                create_vacancy(atomobj) # Create a vacancy

            # Set the momenta corresponding to desired temperature when running Verlocity Verlet
            if Velocity_Verlet_Check == True:
                MaxwellBoltzmannDistribution(atomobj, Temperature * units.kB)
                Stationary(atomobj) # Set linear momentum to zero
                ZeroRotation(atomobj) # Set angular momentum to zero
            # Interatomic potential
            if (EMT_Check == True) and (openKIM_Check == False) and (Lennard_Jones_Check == False):
                atomobj.calc = EMT()
            elif (EMT_Check == False) and (openKIM_Check == True) and (Lennard_Jones_Check == False):
                #Sets the potential for openKIM. If none is given returns standard Lennard-Jones

                potential = checkKIMpotential(KIM_potential)
                #atomobj.calc = KIM(potential, options={"ase_neigh": True})
                atomobj.set_calculator(OpenKIMcalculator(potential))


            elif (EMT_Check == False) and (openKIM_Check == False) and (Lennard_Jones_Check == True):
                atomobj.calc = LennardJones(list(dict.fromkeys(atomobj.get_atomic_numbers())), LJ_epsilon, LJ_sigma, rCut=LJ_cutoff, modified=True)

            else:
                raise Exception("Only one of EMT, OpenKim and Lennard_jones can be true!")

            if (Velocity_Verlet_Check == True) and (Langevin_Check == False):

                # We want to run MD with constant energy using the VelocityVerlet algorithm.
                dyn = VelocityVerlet(atomobj, Time_step*units.fs)
                
            elif (Velocity_Verlet_Check == False) and (Langevin_Check == True):
                dyn = Langevin(atomobj, Time_step*units.fs, units.kB*Temperature, Langevin_friction)
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
            traj.close()

            eq_index = calc.eq_test(atomobj, traj)
            
            if eq_index != 0:#If-statement that checks if we ever reached equilibrium. Returns a message if the traj-file is empty, otherwise does calculations.
                meansSquareDisplacement = calc.MSD_calc(atomobj, traj, -1, eq_index)
                print("MSD = ", meansSquareDisplacement, "[Å²]")
                
                selfDiffusionCoffecient = calc.Self_diffuse(meansSquareDisplacement, (len(traj) - eq_index), Interval, Time_step)
                print("D = ", selfDiffusionCoffecient, "[Å²/fs]")
                
                lindemann = calc.Lindemann(traj, meansSquareDisplacement)
                print("Lindeman Criterion = ", lindemann)
                
                if Velocity_Verlet_Check == True:    
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

                if Run_Optimade == True:
                    translate_to_optimade(atomobj, meansSquareDisplacement, selfDiffusionCoffecient, lindemann , specificHeatCapacity, 
                                            internalTemperature, cohesiveEnergy, internalPressure, B_GPa)

                    concatenateOptimadeDataFiles(Optimade_name) ### Move me to supercomputer script later!!! 
                
                #Writes a .csv-file with time evolution of the mean square displacement
                calc.write_time_evolution_to_csv(atomobj, "Visualization/properties.csv", traj, eq_index, Interval, Time_step)

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

                if Velocity_Verlet_Check == True:
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
                                                    latticeConstant_c,
                                                    Steps*Time_step)


            else:
                print("System never reached equilibrium. No calculations are possible.")

            #Moves the trajectory file to another folder after it has been used
            shutil.move(trajFileName, "Traj/" + trajFileName)
        except:
            print("dyn.run went wrong")
            continue
    atoms.clear()
