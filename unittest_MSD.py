import sys, unittest
from MSD import MSD_calc

from asap3 import EMT
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
dyn = VelocityVerlet(atoms, 5 * units.fs, "atoms.traj")  # 5 fs time step.

class MSDTests(unittest.TestCase):

      def test_MSD_calc(self):
      	  self.assertTrue(True)

if __name__ == '__main__':
   tests = [unittest.TestLoader().loadTestsFromTestCase(MSDTests)]
   testsuite = unittest.TestSuite(tests)
   result = unittest.TextTestRunner(verbosity=0).run(testsuite)
   sys.exit(not result.wasSuccessful())
 
