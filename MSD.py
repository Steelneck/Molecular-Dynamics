"""Calculator for time averaged Mean Square Displacement"""
"""Unclear if MSD should be calculated from t after equilibrium or from equilibrium to t=0"""

from ase.lattice.cubic import FaceCenteredCubic
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.md.verlet import VelocityVerlet
from ase import units
import numpy as np
import tkinter

# Use Asap for a huge performance increase if it is installed
from asap3 import EMT
from asap3 import Trajectory
from asap3.io.trajectory import *

size = 10


    
# Set up a crystal
atoms = FaceCenteredCubic(directions=[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
                          symbol="Cu",
                          size=(size, size, size),
                          pbc=True)

# Describe the interatomic interactions with the Effective Medium Theory
atoms.calc = EMT()

# Set the momenta corresponding to T=300K
MaxwellBoltzmannDistribution(atoms, 300 * units.kB)

# We want to run MD with constant energy using the VelocityVerlet algorithm.
dyn = VelocityVerlet(atoms, 5 * units.fs, "atoms.traj")  # 5 fs time step.


def printenergy(a=atoms):  # store a reference to atoms in the definition.
    """Function to print the potential, kinetic and total energy."""
    epot = a.get_potential_energy() / len(a)
    ekin = a.get_kinetic_energy() / len(a)
    print('Energy per atom: Epot = %.3feV  Ekin = %.3feV (T=%3.0fK)  '
          'Etot = %.3feV' % (epot, ekin, ekin / (1.5 * units.kB), epot + ekin))

def MSD_calc(a): 
    """Function to print the length of the vector to an atom."""
    traj_MSD = Trajectory("atoms.traj")
    pos_eq = traj_MSD[100].get_positions()
    pos_t = traj_MSD[1].get_positions()
    diff = pos_t - pos_eq
    diff_sq = np.absolute(diff)**2
    MSD = np.sum(diff_sq)/len(a)
    print(MSD)

# Now run the dynamics
#traj = Trajectory("atoms.traj", "w", atoms)

#dyn.attach(MSD_calc(atoms), interval=10)
dyn.run(200)
MSD_calc(atoms)
