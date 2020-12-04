""" Unit test for init functions class """

import sys, unittest
import numpy as np
from init_functions import set_lattice_const
from init_functions import set_lattice
from init_functions import insert_impurity
from init_functions import create_vacancy
from ase.lattice.cubic import FaceCenteredCubic
from asap3 import EMT
from ase.build import bulk
from ase import *

atoms = set_lattice(FaceCenteredCubic,
                    0,
                    Directions = [[1,0,0],[0,1,0],[0,0,1]],
                    Miller = [None,None,None],
                    Size_X = 2,
                    Size_Y = 2,
                    Size_Z = 2,
                    Symbol = 'Cu',
                    Pbc = True)

al_cube = bulk('Al', 'fcc', a=3.6)
cu_cube = bulk('Cu', 'fcc', a=3.6)
ag_cube = bulk('Ag', 'fcc', a=3.6)
super_al = al_cube*(2,4,4)
super_cu = cu_cube*(2,4,4)
super_ag = ag_cube*(2,4,4)

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

        def test_set_lattice(self):
            self.assertEqual(atoms.get_chemical_symbols()[0], 'Cu')
            self.assertEqual(len(super_al), len(atoms)) # 32 in each

        def test_atoms_size(self):
            copy_atoms = atoms
            insert_impurity(atoms, 'Au')
            self.assertEqual(len(atoms), len(copy_atoms))

        def test_cube_size(self):
            copy_cube = super_al
            insert_impurity(super_al, 'Ag')
            self.assertEqual(len(copy_cube), len(super_al))

        def test_last_atom(self):
            last_ele_before = atoms.get_chemical_symbols()
            insert_impurity(atoms, 'Ag')
            self.assertIsNot(last_ele_before[-1], atoms.get_chemical_symbols()[-1])

        def test_create_vacancy(self):
            before = len(super_cu)
            create_vacancy(super_cu)
            after = len(super_cu)
            self.assertIsNot(before, after)

        def test_impurity_position(self):
            before = super_ag.get_positions()[0]
            insert_impurity(super_ag, 'Ag')
            after = super_ag.get_positions()[-1]        # Insert_impurity appends the atom at the end of the list
            self.assertEqual(before.all(), after.all()) # Checks if the removed position was the same as insert



        

if __name__ == '__main__':
    tests = [unittest.TestLoader().loadTestsFromTestCase(TestInit)]
    testsuite = unittest.TestSuite(tests)
    result = unittest.TextTestRunner(verbosity=0).run(testsuite)
    sys.exit(not result.wasSuccessful())