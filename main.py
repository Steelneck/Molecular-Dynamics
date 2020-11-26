"""Demonstrates molecular dynamics with constant energy."""

from Init.init_values import *
#import tkinter
import os
import Calculations.calculations as calc
from asap3 import Trajectory

def main():  
    # Initiate the crystal based on the chosen variables
    # This will eventually become "Initiate the system" => system depends on user's choice
    atoms = init()
    
    # We want to run MD with constant energy using the VelocityVerlet algorithm.
    dyn = VelocityVerlet(atoms, 5*units.fs)  # 5 fs time step.

    traj = Trajectory("atoms.traj", "w", atoms)

    dyn.attach(traj.write, interval=10)
    dyn.run(200)

    traj = Trajectory("atoms.traj")
    traj_eq = Trajectory("atoms_eq.traj", "w", atoms)
    
    calc.eq_traj(atoms, traj, traj_eq, Size_X * Size_Y * Size_Z)#Creates new .traj-file containing trajectory post equilibrium.
    if os.path.getsize("atoms_eq.traj") != 0: #If-statement that checks if we ever reached equilibrium. Returns a message if the traj-file is empty, otherwise does calculations.
        traj_eq = Trajectory("atoms_eq.traj")
        MSD = calc.MSD_calc(atoms, traj_eq, 10)
        D = calc.Self_diffuse(traj_eq, MSD, 10)
        L = calc.Lindemann(traj_eq, MSD, 10)
        SHC = calc.Specific_Heat(atoms, traj_eq)
        internalPressure = calc.calc_internal_pressure(atoms, traj_eq, Size_X * Size_Y * Size_Z)
    else:
        print("Something went wront, your system never reached equilibrium. No calculations are possible.")

if __name__ == "__main__":
    main()
