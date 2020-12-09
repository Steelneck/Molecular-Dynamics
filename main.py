"""Demonstrates molecular dynamics with constant energy."""

from Init.init_values import *
from tkinter import *
import os
import Calculations.calculations as calc
from asap3 import Trajectory
from ase import *
import shutil

""" Ignore the parallel rank deprecation warning """

def warn(*args, **kwargs):
    pass
import warnings
warnings.warn = warn

""" Main """

def main():  
    # Initiate the crystal based on the chosen variables
    # This will eventually become "Initiate the system" => system depends on user's choice
    atoms = init()

    # We want to run MD with constant energy using the VelocityVerlet algorithm.
    dyn = VelocityVerlet(atoms, 5*units.fs)  # 5 fs time step.

    traj = Trajectory("atoms.traj", "w", atoms)

    dyn.attach(traj.write, interval)
    dyn.run(steps)
    
    traj = Trajectory("atoms.traj")
    traj_eq = Trajectory("atoms_eq.traj", "w", atoms)

    # I believe that when assigning Calculator to atoms.calc the Calculator object is modified like it is a pointer object/referenced. Thus passing Calculator does not work. 
    # Please feel free to switch EMT() argument to Calculator and see the error message. print(atoms.calc, Calculator) then print(EMT()) proves that Calculator is modified, like a pointer object.
    latticeConstant_a = calc.calc_lattice_constant_fcc_cubic(Symbol, EMT())         
    
    calc.eq_traj(atoms, traj, traj_eq, Size_X * Size_Y * Size_Z)#Creates new .traj-file containing trajectory post equilibrium.
    if os.path.getsize("atoms_eq.traj") != 0: #If-statement that checks if we ever reached equilibrium. Returns a message if the traj-file is empty, otherwise does calculations.
        traj_eq = Trajectory("atoms_eq.traj")
        calc.write_atom_properties(atoms, "Visualization/properties.csv", traj_eq)
        #If-statement that checks if we ever reached equilibrium.
        MSD = calc.MSD_calc(atoms, traj_eq, 1)
        print("MSD = ", MSD, "[Å²]")
        D = calc.Self_diffuse(MSD, len(traj_eq))
        print("D = ", D, "[Å²/fs]")
        L = calc.Lindemann(traj_eq, MSD)
        SHC = calc.Specific_Heat(atoms, traj_eq)

        # Internal temperature of the system
        internalTemperature = calc.internal_temperature(atoms, traj_eq, timeStepIndex)
        
        # Debye temperature of the system
        #debyeTemperature = calc.debye_temperature(atoms, "atoms_eq.traj", 10, EMT())

        internalPressure = calc.calc_internal_pressure(atoms, traj_eq, Size_X * Size_Y * Size_Z)
    else:
        print("System never reached equilibrium. No calculations are possible.")

if __name__ == "__main__":
    main()
