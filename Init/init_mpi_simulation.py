# Initiate the crystal based on the chosen variables
# This will eventually become "Initiate the system" => system depends on user's choice
import os
import shutil
from .init_values import *

# 6 of the 7 lattice systems (rhombohedral is not available)
from ase.lattice.cubic import *
from ase.lattice.tetragonal import *
from ase.lattice.orthorhombic import *
from ase.lattice.monoclinic import *
from ase.lattice.triclinic import *
from ase.lattice.hexagonal import *

from tkinter import *
import Calculations.calculations as calc
from asap3 import Trajectory
from ase.gui import *
from Init.init_functions import set_lattice_const, set_lattice

from mpi4py import MPI
import os
import numpy as np


def run_config(input_config):

    Temperature = input_config[6]
    Interval = input_config[8]
    Steps = input_config[7]
    Size_X = input_config[9]
    Size_Y = input_config[10]
    Size_Z = input_config[11]

    Lattice_Const = set_lattice_const(lc_a = input_config[17],
                                    lc_b = input_config[18],
                                    lc_c = input_config[19],
                                    lc_alpha = input_config[20],
                                    lc_beta = input_config[21],
                                    lc_gamma = input_config[22])

    # Set up a crystals
    
    atoms = set_lattice(getattr(ase.lattice.cubic, input_config[14]),
                        Lattice_Const,
                        Directions = input_config[15],
                        Miller = input_config[16],
                        Size_X = input_config[9],
                        Size_Y = input_config[10],
                        Size_Z = input_config[11],
                        Symbol = input_config[3],
                        Pbc = input_config[13])


    atoms.calc = getattr(asap3, input_config[0])() # = EMT()
    
        # Set the momenta corresponding to desired temperature when running Verlocity Verlet
    if input_config[1] == 'VelocityVerlet':
        MaxwellBoltzmannDistribution(atoms, Temperature * units.kB)

    Stationary(atoms) # Set linear momentum to zero

    dyn = getattr(asap3.md.verlet, input_config[1])(atoms, 5*units.fs)

    #Creates a unique name for every simulation run 
    trajFileName = atoms.get_chemical_formula() + '.traj'
    traj = Trajectory(trajFileName, "w", atoms)
    dyn.attach(traj.write, Interval)
    dyn.run(Steps)
    
    traj = Trajectory(trajFileName)

    #latticeConstant_a = calc.calc_lattice_constant_fcc_cubic(Symbol, EMT())
    #print("Lattice constant a:", latticeConstant_a)
    
    eq_index = calc.eq_traj(atoms, traj, Size_X * Size_Y * Size_Z)
    if eq_index != 0:
        
        MSD = calc.MSD_calc(atoms, traj, eq_index)
        print("MSD = ", MSD, "[Å²]")
        
        D = calc.Self_diffuse(MSD, (len(traj) - eq_index))
        print("D = ", D, "[Å²/fs]")
        
        L = calc.Lindemann(traj, MSD)
        
        SHC = calc.Specific_Heat(atoms, traj, eq_index)
        print("C_p = ", SHC, "[J/K*Kg]")
        
        internalTemperature = calc.internal_temperature(atoms, traj, eq_index)
        print("Internal temperature:", internalTemperature, "[K]")
        
        cohesiveEnergy = calc.cohesive_energy(atoms, traj, eq_index)
        print("Cohesive energy:", cohesiveEnergy, "[eV/atom]")
        
        internalPressure = calc.calc_internal_pressure(atoms, traj, eq_index, Size_X * Size_Y * Size_Z)
        print("Internal Pressure:", internalPressure, "[eV / Å^3]")
        
        # e0, v0, B_GPa = calc.calc_bulk_modulus(atoms)
        # print('Bulk Modulus:', B_GPa, '[GPa]', '|', 'Minimum energy E =', e0, '[eV], at volume V =', v0, '[Å^3].')
    else:
        print("System never reached equilibrium. No calculations are possible.")

    #Moves the trajectory file to another folder after it has been used
    #shutil.move(trajFileName, "Traj/" + trajFileName)


# def simulation_mpi(input_list):


#     # Här fixar jag query

    # for input_config in input_list:
    #     run_config(input_config)

