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

    traj_MSD = Trajectory("atoms.traj")
    n = 1
    atoms_eq = []
    while n < len(traj_MSD):
        ediff = abs((traj_MSD[n].get_potential_energy() / len(atoms)) -
               (traj_MSD[n].get_kinetic_energy() / len(atoms))) #epot-ekin
        print(atoms.get_potential_energy() / len(atoms))
        print(atoms.get_kinetic_energy() / len(atoms))
        if ediff < 0.005:
            atoms_eq.append(traj_MSD[n])
        n += 1
    SH = calc.Specific_Heat(atoms, atoms_eq)
    #MSD = calc.MSD_calc(atoms, 10, atoms_eq)
    #D = calc.Self_diffuse(MSD, 10, atoms_eq)
    #L = calc.Lindemann(atoms, MSD)

if __name__ == "__main__":
    main()
