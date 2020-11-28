"""Demonstrates molecular dynamics with constant energy."""

from Init.init2 import *
import os
from tkinter import *
import Calculations.calculations as calc
from asap3 import Trajectory

def main():

    # Initiate the crystal based on the chosen variables
    # This will eventually become "Initiate the system" => system depends on user's choice
    
    atoms = init_MP()
    
    
    # We want to run MD with constant energy using the VelocityVerlet algorithm.
    for atomobj in atoms:
        
        dyn = VelocityVerlet(atomobj, 5 * units.fs)  # 5 fs time step.

        traj = Trajectory("atoms.traj", "w", atomobj)
        
        dyn.attach(traj.write, interval=10)
        dyn.run(100)
        
        calc.eq_traj(atomobj) #Creates new .traj-file containing trajectory post equilibrium.
        if os.path.getsize("atoms_eq.traj") != 0:
            traj_eq = Trajectory("atoms_eq.traj")
            #If-statement that checks if we ever reached equilibrium.
            MSD = calc.MSD_calc(atomobj, 10)
            D = calc.Self_diffuse(MSD, 10)
            L = calc.Lindemann(atomobj, MSD)
            SHC = calc.Specific_Heat(atomobj)
        else:
            print("Something went wront, your system never reached equilibrium. No calculations are possible.")

if __name__ == "__main__":
    main()
