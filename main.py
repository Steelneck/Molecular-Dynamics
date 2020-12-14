"""Demonstrates molecular dynamics with constant energy."""

from Init.init_values import *
from tkinter import *
import os
import Calculations.calculations as calc
from asap3 import Trajectory
from ase import *
#from ase.geometry import crystal_structure_from_cell, Cell
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
        print("Mean Square Displacement = ", MSD, "[Å²]")
        D = calc.Self_diffuse(MSD, len(traj_eq))
        print("Self Diffusion coefficient = ", D, "[Å²/fs]")
        L = calc.Lindemann(traj_eq, MSD)
        print("L = ", L)
        SHC = calc.Specific_Heat(atoms, traj_eq)

        # Internal temperature of the system
        Temp = calc.internal_temperature(atoms, traj_eq)
        print("Internal temperature: T =", Temp, "[K]")

        # Cohesive energy of the system
        Ecoh = calc.cohesive_energy(atoms, traj_eq)
        print("Cohesive energy: Eᴄᴏʜ =", Ecoh, "[eV/atom]")

        # Debye temperature of the system
        Debye = calc.debye_temperature(traj_eq, MSD)
        print("Debye temperature: Θ =", Debye, "[K]")

        internalPressure = calc.calc_internal_pressure(atoms, traj_eq, Size_X * Size_Y * Size_Z)
        e0, v0, B_GPa = calc.calc_bulk_modulus(atoms)
        
        calc.write_simulation_values(Symbol,
                                     atoms.get_cell().get_bravais_lattice(),
                                     Size_X * Size_Y * Size_Z,
                                     latticeConstant_a,
                                     Temperature,
                                     MSD,
                                     D,
                                     L,
                                     SHC,
                                     internalTemperature,
                                     internalPressure,
                                     cohesiveEnergy)
    else:
        print("System never reached equilibrium. No calculations are possible.")

if __name__ == "__main__":
    main()
