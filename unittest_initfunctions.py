""" Unit test for init functions class """

import sys, unittest
from init_functions import set_lattice_const
from init_functions import set_lattice
from ase.lattice.cubic import FaceCenteredCubic
from asap3 import EMT

" Define the test class with test functions "

class TestInit(unittest.TestCase):
        
        def test_lc_none(self):
            lc_constant_zero = set_lattice_const(0, 0, 0, 0, 0, 0)
            self.assertEqual(lc_constant_zero, None)
            
        def test_lc_zero(self):
            lc_constant_zero = set_lattice_const(0, 0, 0, 0, 0, 0)
            self.assertTrue(lc_constant_zero != 0)
            
        def test_lc_a(self):
            lc_constant_a = set_lattice_const(3.6, 0, 0, 0, 0, 0)
            self.assertTrue(isinstance(lc_constant_a, float))
            
        def test_lc_dict(self):
            lc_constant_dict = set_lattice_const(3.6, 0, 3.8, 0, 0, 74)
            self.assertTrue(isinstance(lc_constant_dict, dict))
            
        def test_lc_values_exist(self):
            lc_constant = set_lattice_const(3.6, 3.8, 2, 23, 45, 74)
            self.assertTrue('a' in lc_constant)
            self.assertTrue('b' in lc_constant)
            self.assertTrue('c' in lc_constant)
            self.assertTrue('alpha' in lc_constant)
            self.assertTrue('beta' in lc_constant)
            self.assertTrue('gamma' in lc_constant)
            
        def test_lc_values(self):
            lc_constant = set_lattice_const(3.6, 0, 2, 0, 0, 74)
            self.assertEqual(lc_constant['a'], 3.6)
            self.assertEqual(lc_constant['c'], 2)
            self.assertEqual(lc_constant['gamma'], 74)

if __name__ == '__main__':
    tests = [unittest.TestLoader().loadTestsFromTestCase(TestInit)]
    testsuite = unittest.TestSuite(tests)
    result = unittest.TextTestRunner(verbosity=0).run(testsuite)
    sys.exit(not result.wasSuccessful())