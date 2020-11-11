import sys, unittest
from ase.lattice.cubic import FaceCenteredCubic
from calculations import Specific_Heat
import numpy
""" 
    # Example From hands-on-3, delete later
from ase.lattice.cubic import FaceCenteredCubic
from ase.calculators.emt import EMT
"""

atoms = FaceCenteredCubic(directions=[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
                                  symbol="Cu",
                                  size=(3, 3, 3),
                                  pbc=True)

not_bravice_lattice = 1


class PropertyCalculationTests(unittest.TestCase):
 

    def test_specific_heat(self):
        self.assertIsInstance(Specific_Heat(atoms), numpy.float64)
        
    def test_specific_heat_not_bravice_lattice(self):
        self.assertIsNone(Specific_Heat(not_bravice_lattice))
        


if __name__ == '__main__':
    tests = [unittest.TestLoader().loadTestsFromTestCase(PropertyCalculationTests)]
    testsuite = unittest.TestSuite(tests)
    result = unittest.TextTestRunner(verbosity=0).run(testsuite)
    sys.exit(not result.wasSuccessful())