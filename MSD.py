"""Calculator for time averaged Mean Square Displacement"""
"""Unclear if MSD should be calculated from t after equilibrium or from equilibrium to t=0"""

from ase.lattice.cubic import FaceCenteredCubic
from ase.lattice.cubic import BodyCenteredCubic
from ase.lattice.cubic import SimpleCubic
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.md.verlet import VelocityVerlet
from ase.lattice.bravais import Lattice
from ase.atoms import Atoms
from ase.md.langevin import Langevin
from ase import units

import numpy as np
import tkinter

# Use Asap for a huge performance increase if it is installed
from asap3 import EMT, Trajectory, FullNeighborList
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
dyn = VelocityVerlet(atoms, 5 * units.fs)  # 5 fs time step.
# dyn = Langevin(atoms, 5 * units.fs, units.kB * 300, 0.002)

def MSD_calc(a, t): 
    """Function to calculate and print the time average of the mean square displacement (MSD) at time t. Also calculates and prints the self-diffusion coefficient(D)."""
    traj_MSD = Trajectory("atoms.traj")  #Fetch atom trajectories from file.
    time = len(traj_MSD)-t
    pos_eq = traj_MSD[100].get_positions() #position of all atoms at time after equilibrium has been reached.
    pos_t = traj_MSD[time].get_positions() #position of all atoms at time t
    diff = pos_t - pos_eq #displacement of all atoms from equilibrium to time t as a vector
    diff_sq = np.absolute(diff)**2 #Square of the displacement
    MSD = np.sum(diff_sq)/len(a) #Time averaged mean square displacement.
    D = MSD/(6*time) #How to connect mean squre displacement to self-diffusion coefficient.
    open("atoms.traj", "w").close()
    #print("MSD = ", MSD)
    #print("D = ", D)
    return(MSD, D)

def Lindemann(a, MSD):
     nblist = FullNeighborList(6, atoms).get_neighbors(1, -1) #Returns a list containing information about neighbors to atom 1 at a maximum distance of 6Å. Third element of nblist is a list of all squared norms of distances between atom 1 and its neighbors.
     d = np.sqrt(np.amin(nblist[2])) #d is the distance to the nearest neighbor. nblist[2] since the third element of nblist is a list of the square of the norm of all distances.
     L = MSD/d #Lindemann criterion. Expect melting when L>0.1
     if L > 0.1:
         #print("MELTING!")
         return True
     else:
         #print("Not melting.")
         return False




# Now run the dynamics

traj = Trajectory("atoms.traj", "w", atoms)
dyn.attach(traj.write, interval=1)
dyn.run(200)

MSD = MSD_calc(atoms, 50)
L = Lindemann(atoms, MSD[1])
