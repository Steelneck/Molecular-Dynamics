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
    dyn = VelocityVerlet(atoms, 5 * units.fs)  # 5 fs time step.

    traj = Trajectory("atoms.traj", "w", atoms)
    
    dyn.attach(traj.write, interval=10)
    dyn.run(200)
    
    calc.eq_traj(atoms) #Creates new .traj-file containing trajectory post equilibrium.
    if os.path.getsize("atoms_eq.traj") != 0:
        traj_eq = Trajectory("atoms_eq.traj")
        #If-statement that checks if we ever reached equilibrium.
        MSD = calc.MSD_calc(atoms, 10)
        D = calc.Self_diffuse(MSD, 10)
        L = calc.Lindemann(atoms, MSD)
        SHC = calc.Specific_Heat(atoms)
        internalPressure = calc.calc_internal_pressure(atoms, traj_eq)
    else:
        print("Something went wront, your system never reached equilibrium. No calculations are possible.")

if __name__ == "__main__":
    main()
