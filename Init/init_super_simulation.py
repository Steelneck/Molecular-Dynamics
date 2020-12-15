# Initiate the crystal based on the chosen variables
# This will eventually become "Initiate the system" => system depends on user's choice

from Init.init_values import *
from tkinter import *
import os
from tkinter import *
import Calculations.calculations as calc
from asap3 import Trajectory
from ase.gui import *

def super_simulation(Calculator,Symbol):
    """ Choose which init function to run. 
            init() for ASE configuration.
            init_MP() for materials project configuration.
        See init_values for configuration settings.
    """
    atoms = init(Calculator,Symbol)
    
    for atomobj in atoms:
        # We want to run MD with constant energy using the VelocityVerlet algorithm.
        dyn = VelocityVerlet(atomobj, 5*units.fs)  # 5 fs time step.

        fileName = Symbol + ".traj"
        traj = Trajectory(fileName, "w", atomobj)

        dyn.attach(traj.write, interval)
        dyn.run(steps)
        
        traj = Trajectory(fileName)

        eq_index = calc.eq_traj(atomobj, traj, Size_X * Size_Y * Size_Z)
        print(eq_index)
        print("Calc of", Symbol, "gave:")
        
        #If-statement that checks if we ever reached equilibrium.
        if eq_index != 0:
            MSD = calc.MSD_calc(atomobj, traj, eq_index)
            D = calc.Self_diffuse(traj, eq_index, MSD)
            L = calc.Lindemann(traj, MSD)
            SHC = calc.Specific_Heat(atomobj, traj, eq_index)

            # Internal temperature of the system
            internalTemperature = calc.internal_temperature(atomobj, traj)

        else:
            print("System never reached equilibrium. No calculations are possible.")
    atoms.clear()
