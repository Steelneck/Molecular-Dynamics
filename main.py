"""Demonstrates molecular dynamics with constant energy."""

from Init.init_values import *
#import tkinter
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

    
    calc.eq_traj(atoms)
    MSD = calc.MSD_calc(atoms, 10)
    D = calc.Self_diffuse(MSD, 10)
    L = calc.Lindemann(atoms, MSD)

if __name__ == "__main__":
    main()
