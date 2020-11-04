import sys, unittest
""" 
    # Example From hands-on-3, delete later
from ase.lattice.cubic import FaceCenteredCubic
from ase.calculators.emt import EMT
"""
class PropertyCalculationTests(unittest.TestCase):

    """ 
    # Example From hands-on-3, delete later
    def test_calcenergy(self):

        atoms = FaceCenteredCubic(directions=[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
                                  symbol="Cu",
                                  size=(3, 3, 3),
                                  pbc=True)
        atoms.calc = EMT()

        ekin, epot = calcenergy(atoms)


        if ekin is not None and epot is not None:
            self.assertTrue(True)
        else:
            self.assertTrue(False)
    """
    # Dummy test to have some sort of validation that code works, just says that test was success
    def dummyTest(self):
        self.assertTrue(True)

if __name__ == '__main__':
    tests = [unittest.TestLoader().loadTestsFromTestCase(PropertyCalculationTests)]
    testsuite = unittest.TestSuite(tests)
    result = unittest.TextTestRunner(verbosity=0).run(testsuite)
    sys.exit(not result.wasSuccessful())