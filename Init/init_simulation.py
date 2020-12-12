# Initiate the crystal based on the chosen variables
# This will eventually become "Initiate the system" => system depends on user's choice
import os
import shutil
from Init.init_values import *
from tkinter import *
from tkinter import *
import Calculations.calculations as calc
from asap3 import Trajectory
from ase.gui import *


def simulation(EMT_Check,openKIM_Check,KIM_potential,ASE, Materials_project,Symbol,Critera_list, 
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
                atoms = init(EMT_Check, openKIM_Check, KIM_potential,Symbol,
                            Vacancy, Impurity, Impurity_ele, Temperature,
                            Size_X,Size_Y,Size_Z,PBC,Directions,Miller,
                            lc_a,lc_b,lc_c,lc_alpha,lc_beta,lc_gamma)

        else:
            atoms = init(EMT_Check, openKIM_Check, KIM_potential,Symbol,
                            Vacancy, Impurity, Impurity_ele_list, Temperature,
                            Size_X,Size_Y,Size_Z,PBC,Directions,Miller,
                            lc_a,lc_b,lc_c,lc_alpha,lc_beta,lc_gamma)
    elif (Materials_project == True) and (ASE == False):
        if Impurity == True:
            for Impurity_ele in Impurity_ele_list:
                atoms = init_MP(EMT_Check,openKIM_Check,KIM_potential,Critera_list,
                                Vacancy, Impurity, Impurity_ele, Temperature,
                                Size_X,Size_Y,Size_Z,API_Key,PBC)

        else:
            atoms = init_MP(EMT_Check,openKIM_Check,KIM_potential,Critera_list,
                                Vacancy, Impurity, Impurity_ele_list, Temperature,
                                Size_X,Size_Y,Size_Z,API_Key,PBC)
    else:
        raise Exception("ASE=Materials_Materials. Both cannot be true/false at the same time!")
    
    for atomobj in atoms:

        # We want to run MD with constant energy using the VelocityVerlet algorithm.
        dyn = VelocityVerlet(atomobj, 5*units.fs)  # 5 fs time step.

        #Creates a unique name for every simulation run 
        trajFileName = atomobj.get_chemical_formula() + '.traj'
        trajFileName_eq = atomobj.get_chemical_formula() + '_eq.traj'
        traj = Trajectory(trajFileName, "w", atomobj)
        
        dyn.attach(traj.write, Interval)
        dyn.run(Steps)
        
        traj = Trajectory(trajFileName)
        traj_eq = Trajectory(trajFileName_eq, "w", atomobj)

        latticeConstant_a = calc.calc_lattice_constant_fcc_cubic(Symbol, EMT())
        
        calc.eq_traj(atomobj, traj, traj_eq, Size_X * Size_Y * Size_Z)#Creates new .traj-file containing trajectory post equilibrium.
        if os.path.getsize(trajFileName_eq) != 0: #If-statement that checks if we ever reached equilibrium. Returns a message if the traj-file is empty, otherwise does calculations.
            traj_eq = Trajectory(trajFileName_eq)
            calc.write_atom_properties(atoms, "Visualization/properties.csv", traj_eq)
            
            MSD = calc.MSD_calc(atomobj, traj_eq, 1)
            print("MSD = ", MSD, "[Å²]")
            
            D = calc.Self_diffuse(MSD, len(traj_eq))
            print("D = ", D, "[Å²/fs]")
            
            L = calc.Lindemann(traj_eq, MSD)
            
            SHC = calc.Specific_Heat(atomobj, traj_eq)
            print("C_p = ", SHC, "[J/K*Kg]")
            
            internalTemperature = calc.internal_temperature(atomobj, traj_eq)
            print("Internal temperature:", internalTemperature, "[K]")
            
            cohesiveEnergy = calc.cohesive_energy(atomobj, traj_eq)
            print("Cohesive energy:", cohesiveEnergy, "[eV/atom]")
            
            internalPressure = calc.calc_internal_pressure(atomobj, traj_eq, Size_X * Size_Y * Size_Z)
            print("Internal Pressure:", internalPressure, "[eV / Å^3]")
            
            e0, v0, B_GPa = calc.calc_bulk_modulus(atomobj)
            print('Bulk Modulus:', B_GPa, '[GPa]', '|', 'Minimum energy E =', e0, '[eV], at volume V =', v0, '[Å^3].')

            #Moves the trajectory file to another folder after it has been used
            shutil.move(trajFileName, "Traj/" + trajFileName)
            shutil.move(trajFileName_eq, "Traj/" + trajFileName_eq)

        else:
            print("System never reached equilibrium. No calculations are possible.")
    atoms.clear()
