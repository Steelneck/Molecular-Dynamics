import sys, unittest
from MSD import MSD_calc
from MSD import Lindemann
from MSD import Self_diffuse

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

from asap3 import EMT, Trajectory, FullNeighborList
from asap3 import Trajectory
from asap3.io.trajectory import *

size = 10

# Set up a crystal to use as test object

atoms = FaceCenteredCubic(directions=[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
                          symbol="Cu",
                          size=(size, size, size),
                          pbc=True)
			  
atoms.calc = EMT()
MaxwellBoltzmannDistribution(atoms, 300 * units.kB)
dyn = VelocityVerlet(atoms, 5 * units.fs)  # 5 fs time step.

traj = Trajectory("atoms.traj", "w", atoms)
dyn.attach(traj.write, interval=1)
dyn.run(200)

traj_MSD = Trajectory("atoms.traj")
n = 1
atoms_eq = []
while n < len(traj_MSD):
      ediff = abs((traj_MSD[n].get_potential_energy() / len(atoms)) -
               (traj_MSD[n].get_kinetic_energy() / len(atoms))) #epot-ekin
      if ediff < 0.005:
            atoms_eq.append(traj_MSD[n])
      n += 1

MSD = MSD_calc(atoms, 19, atoms_eq)
D = Self_diffuse(MSD, 19, atoms_eq)
L = Lindemann(atoms, MSD)

open("atoms.traj", "w").close()

class MSDTests(unittest.TestCase):
      def test_MSD_calc(self):
      	    self.assertTrue(isinstance(MSD, float))
      def test_self_diffuse(self):
            self.assertTrue(isinstance(D, float))
      def test_Lindemann(self):
            self.assertTrue(isinstance(L, int))

if __name__ == '__main__':
   tests = [unittest.TestLoader().loadTestsFromTestCase(MSDTests)]
   testsuite = unittest.TestSuite(tests)
   result = unittest.TextTestRunner(verbosity=0).run(testsuite)
   sys.exit(not result.wasSuccessful())
 
