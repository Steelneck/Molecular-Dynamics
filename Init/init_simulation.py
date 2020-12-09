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

def simulation(Name, EMT_Check,openKIM_Check,KIM_potential,ASE, Materials_project,
                        Symbol, Temperature, Steps, Interval,
                        Size_X, Size_Y, Size_Z,API_Key,PBC,Directions,Miller,
                        lc_a,lc_b,lc_c,lc_alpha,lc_beta,lc_gamma):
    """ Choose which init function to run. 
            init() for ASE configuration.
            init_MP() for materials project configuration.
        See init_values for configuration settings.
    """
    if (ASE == True) and (Materials_project == False) :
        atoms = init(EMT_Check, openKIM_Check, KIM_potential,Symbol,Temperature,Size_X,Size_Y,Size_Z,PBC,Directions,Miller,
                    lc_a,lc_b,lc_c,lc_alpha,lc_beta,lc_gamma)
    elif (Materials_project == True) and (ASE == False):
        atoms = init_MP(EMT_Check,openKIM_Check,KIM_potential,Symbol,Temperature,Size_X,Size_Y,Size_Z,API_Key,PBC)
    else:
        raise Exception("ASE=Materials_Materials. Both cannot be true/false at the same time!")
    
    for atomobj in atoms:
        timeStepIndex = timestepindex(Steps, Interval)
        # We want to run MD with constant energy using the VelocityVerlet algorithm.
        dyn = VelocityVerlet(atomobj, 5*units.fs)  # 5 fs time step.
        trajFileName = Name + '.traj'
        trajFileName_eq = Name + '_eq.traj'
        traj = Trajectory(trajFileName, "w", atomobj)
        
        dyn.attach(traj.write, Interval)
        dyn.run(Steps)
        
        traj = Trajectory(trajFileName)
        traj_eq = Trajectory(trajFileName_eq, "w", atomobj)
        
        calc.eq_traj(atomobj, traj, traj_eq, Size_X * Size_Y * Size_Z)#Creates new .traj-file containing trajectory post equilibrium.
        if os.path.getsize(trajFileName_eq) != 0: #If-statement that checks if we ever reached equilibrium. Returns a message if the traj-file is empty, otherwise does calculations.
            traj_eq = Trajectory(trajFileName_eq)
            #If-statement that checks if we ever reached equilibrium.
            MSD = calc.MSD_calc(atomobj, traj_eq, timeStepIndex)
            D = calc.Self_diffuse(traj_eq, MSD)
            L = calc.Lindemann(traj_eq, MSD)
            SHC = calc.Specific_Heat(atomobj, traj_eq)

            # Internal temperature of the system
            internalTemperature = calc.internal_temperature(atomobj, traj_eq, timeStepIndex)

            #Moves the trajectory file to another folder for later use
            shutil.move(trajFileName, "Traj/" + trajFileName)
            shutil.move(trajFileName_eq, "Traj/" + trajFileName_eq)

        else:
            print("System never reached equilibrium. No calculations are possible.")
    atoms.clear()
