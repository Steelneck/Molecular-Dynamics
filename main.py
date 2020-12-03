"""Demonstrates molecular dynamics with constant energy."""

from Init.init_values import *
from tkinter import *
import os
from tkinter import *
import Calculations.calculations as calc
from asap3 import Trajectory
from ase.gui import *

def main():  
    # Initiate the crystal based on the chosen variables
    # This will eventually become "Initiate the system" => system depends on user's choice

    atoms = init_MP()
    
    for atomobj in atoms:
        # We want to run MD with constant energy using the VelocityVerlet algorithm.
        dyn = VelocityVerlet(atomobj, 5*units.fs)  # 5 fs time step.


        traj = Trajectory("atoms.traj", "w", atomobj)

        dyn.attach(traj.write, interval)
        dyn.run(steps)
        
        traj = Trajectory("atoms.traj")
        traj_eq = Trajectory("atoms_eq.traj", "w", atomobj)
        
        calc.eq_traj(atomobj, traj, traj_eq, Size_X * Size_Y * Size_Z)#Creates new .traj-file containing trajectory post equilibrium.
        if os.path.getsize("atoms_eq.traj") != 0: #If-statement that checks if we ever reached equilibrium. Returns a message if the traj-file is empty, otherwise does calculations.
            traj_eq = Trajectory("atoms_eq.traj")
            #If-statement that checks if we ever reached equilibrium.
            MSD = calc.MSD_calc(atomobj, traj_eq, timeStepIndex)
            D = calc.Self_diffuse(traj_eq, MSD)
            L = calc.Lindemann(traj_eq, MSD)
            SHC = calc.Specific_Heat(atomobj, traj_eq)

            # Internal temperature of the system
            internalTemperature = calc.internal_temperature(atomobj, traj_eq, timeStepIndex)

        else:
            print("System never reached equilibrium. No calculations are possible.")

if __name__ == "__main__":
    main()
