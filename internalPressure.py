""" """
""" """

from ase.lattice.cubic import FaceCenteredCubic
#from ase.lattice.cubic import BodyCenteredCubic
#from ase.lattice.cubic import SimpleCubic
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.md.verlet import VelocityVerlet
from ase.lattice.bravais import Lattice
from ase.atoms import Atoms
from ase.md.langevin import Langevin
from ase import units
from ase import calculators

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

def calc_instantaneous_pressure(myAtoms, trajectoryFileName, timeStepIndex):
    """ Calculates instantaneous pressure at time timeStepIndex * deltaT, i.e. P(n*deltaT)."""
    traj = Trajectory(trajectoryFileName)
    equilibriumIndex = timeStepIndex
    atomsAtEquilibrium = traj[equilibriumIndex]             # Retrieve atoms some amount of timesteps in when system has reached equilibrium
    eqPos = atomsAtEquilibrium.get_positions()              # Atom 3D coords.
    eqTemperature = atomsAtEquilibrium.get_temperature()    # Temperature at time t = equilibriumIndex * deltaT.
    eqPotEn = atomsAtEquilibrium.get_potential_energy()     # Get potetntial energy from atoms object at time t = equilibriumIndex * deltaT, should be the same for all. 
    N = len(atomsAtEquilibrium)                             # Number of atoms in object
    eqForces = atoms.calc.get_forces(atomsAtEquilibrium)    # Use calculator object to get the forces on the atoms. 
    
    # Instantaneos pressure 
    posForce = np.sum(eqPos * eqForces, axis=1)             # Dot product of all forces and positions
    instantP = (N * units.kB * eqTemperature) / eqPotEn + np.sum(posForce)

    #print(instantP)
    return(instantP)

def calc_internal_pressure(myAtoms, trajectoryFileName, iterations):
    """ Internal pressure is the MD average of the instantaneous pressures """
    M = iterations                                          # M = (iterations * deltaT) / deltaT

    # MD average
    allInstantPressures = 0
    for n in range(1,iterations):
        allInstantPressures += calc_instantaneous_pressure(myAtoms, trajectoryFileName, n)
    
    internalPressure = allInstantPressures / M              # Internal pressure is the MD average of the instantaneous pressures
    print("Internal Pressure", internalPressure)
    return(internalPressure)


# Now run the dynamics
traj = Trajectory("atoms.traj", "w", atoms) # Generate trajectory file
dyn.attach(traj.write, interval=1)
iterations = 200
dyn.run(iterations)

internalPressure = calc_internal_pressure(atoms, "atoms.traj", iterations)
