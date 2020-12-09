import sys, unittest, os

from ase.lattice.cubic import FaceCenteredCubic

from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.md.verlet import VelocityVerlet
from ase import units

from calculations import Specific_Heat
from calculations import internal_temperature
from calculations import cohesive_energy
from calculations import calc_instantaneous_pressure
from calculations import calc_internal_pressure
from calculations import eq_traj
from calculations import MSD_calc
from calculations import Self_diffuse
from calculations import Lindemann
from calculations import calc_lattice_constant_fcc_cubic
from calculations import write_atom_properties

import numpy
from asap3 import Trajectory, EMT

atoms = FaceCenteredCubic(directions=[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
                                  symbol="Cu",
                                  size=(3, 3, 3),
                                  pbc=True)

not_bravice_lattice = 1

#Attach calculator to the atoms object.
atoms.calc = EMT()

MaxwellBoltzmannDistribution(atoms, 300 * units.kB)

# We want to run MD with constant energy using the VelocityVerlet algorithm.
dyn = VelocityVerlet(atoms, 5 * units.fs)  # 5 fs time step.

# Create trajectory object
trajWriter = Trajectory("test.traj", "w", atoms)

dyn.attach(trajWriter.write, interval=10)
dyn.run(200)

trajObject = Trajectory("test.traj")

class PropertyCalculationTests(unittest.TestCase):
 
    """Unittests for specific heat"""
    
    def test_specific_heat(self):
        self.assertIsInstance(Specific_Heat(atoms, trajObject), numpy.float64)
        
    def test_specific_heat_not_bravice_lattice(self):
        self.assertFalse(Specific_Heat(not_bravice_lattice, trajObject))

    """Unittests for instantaneous and internal pressure"""
    
    def test_instantaneous_pressure_return_type(self):
        self.assertIsInstance(calc_instantaneous_pressure(atoms, trajObject, 1000, 1), float)

    def test_instantaneous_pressure_wrong_input_argument(self):
        instPressure1 = calc_instantaneous_pressure(None, trajObject, 1000, 1)
        instPressure2 = calc_instantaneous_pressure(atoms, None, 1000, 1)
        instPressure3 = calc_instantaneous_pressure(atoms, trajObject, None, 1)
        instPressure4 = calc_instantaneous_pressure(atoms, trajObject, 1000, None)
        
        # Then all should return None
        self.assertIsNone(instPressure1)
        self.assertIsNone(instPressure2)
        self.assertIsNone(instPressure3)
        self.assertIsNone(instPressure4)

    def test_internal_pressure_return_type(self):
        self.assertIsInstance(calc_internal_pressure(atoms, trajObject, 1000), float)

    def test_internal_pressure_wrong_input_argument(self):
        internalPressure1 = calc_internal_pressure(None, trajObject, 1000)
        internalPressure2 = calc_internal_pressure(atoms, None, 1000)
        internalPressure3 = calc_internal_pressure(atoms, trajObject, None)
        
        # Then all should return None
        self.assertIsNone(internalPressure1)
        self.assertIsNone(internalPressure2)
        self.assertIsNone(internalPressure3)

    """Unittest for eq_calc"""

    def test_eq_calc_check_traj(self):
        eq_trajObject = Trajectory("test_eq.traj", "w", atoms)
        eq_traj(atoms, trajObject, eq_trajObject, 3*3*3)
        eq_trajObject = Trajectory("test_eq.traj")
        self.assertIsInstance(eq_trajObject[1].get_potential_energy(), float)

    #eq_traj doesnt use atoms yet, so no point in testing that input
    def test_eq_calc_wrong_input_argument(self):
        eq_trajObject = Trajectory("test_eq.traj", "w", atoms)
        eq_traj1 = eq_traj(atoms, None, eq_trajObject, 3*3*3)
        eq_traj2 = eq_traj(atoms, trajObject, None, 3*3*3)
        eq_traj3 = eq_traj(atoms, trajObject, eq_trajObject, None)

        #All should return None
        self.assertIsNone(eq_traj1)
        self.assertIsNone(eq_traj2)
        self.assertIsNone(eq_traj3)

    """Unittests for calculation of mean square displacement"""

    def test_MSD_calc_return_type(self):
      	self.assertIsInstance(MSD_calc(atoms, trajObject, 10), float)

    #MSD doesnt use the time input yet so no point in testing it    
    def test_MSD_calc_wrong_input_argument(self):   
        MSD1 = MSD_calc(None, trajObject, 10)
        MSD2 = MSD_calc(atoms, None, 10)
        MSD3 = MSD_calc(atoms, trajObject, None)

        #All should return None
        self.assertIsNone(MSD1)
        self.assertIsNone(MSD2)
        self.assertIsNone(MSD3)

    """Unittests for calculation of Self diffusion coefficient"""

    def test_self_diffuse_return_type(self):
        self.assertIsInstance(Self_diffuse(MSD_calc(atoms, trajObject, 10), len(trajObject)), float)

    #Self_diffuse doesnt use the time input yet so no point in testing it 
    def test_Self_diffuse_wrong_input_argument(self):
        D1 = Self_diffuse(None, len(trajObject))
        D2 = Self_diffuse(MSD_calc(atoms, trajObject, 10), None)

        #All should return None
        self.assertIsNone(D1)
        self.assertIsNone(D2)
        
    def test_Lindemann_return_type(self):
        self.assertIsInstance(Lindemann(trajObject, MSD_calc(atoms, trajObject, 10)), int)

    def test_internal_temperature(self):
        self.assertIsInstance(internal_temperature(atoms, trajObject), float)

    def test_internal_temperature_not_negative(self):
        self.assertGreaterEqual(internal_temperature(atoms, trajObject), 0)

    def test_cohesive_energy(self):
        self.assertIsInstance(cohesive_energy(atoms, trajObject), float)

    def test_cohesive_energy_positive(self):
        self.assertGreater(cohesive_energy(atoms, trajObject), 0)

    #Lindemann doesnt use the time input yet so no point in testing it 
    def test_Lindemann_wrong_input_argument(self):
        L1 =Lindemann(None, MSD_calc(atoms, trajObject, 10))
        L2 =Lindemann(trajObject, None)

        #All should return None
        self.assertIsNone(L1)
        self.assertIsNone(L2)

    """Unittests for calculation of Lattice Constant"""
    def test_lattice_constant_wrong_input_argument(self):
        a1 = calc_lattice_constant_fcc_cubic(None, EMT())
        a2 = calc_lattice_constant_fcc_cubic("not an atom", EMT())
        a3 = calc_lattice_constant_fcc_cubic("Cu", None)
        
        self.assertIsNone(a1)
        self.assertIsNone(a2)
        self.assertIsNone(a3)
    
    def test_lattice_constant_output_type(self):
        self.assertIsInstance(calc_lattice_constant_fcc_cubic('Cu', EMT()), float)

    
    def test_csv_writer_wrong_input_argument(self):
        csv1 = write_atom_properties(None, "properties_test.csv", trajObject)
        csv2 = write_atom_properties(atoms, None, trajObject)
        csv3 = write_atom_properties(atoms, "properties.csv", None)

        #All should return none.
        self.assertIsNone(csv1)
        self.assertIsNone(csv2)
        self.assertIsNone(csv3)

        
    def test_csv_writer_check_csv(self):
        #Check that the .csv-file exists and is not empty.
        eq_trajObject = Trajectory("test_eq.traj", "w", atoms)
        eq_traj(atoms, trajObject, eq_trajObject, 3*3*3)
        eq_trajObject = Trajectory("test_eq.traj")
        write_atom_properties(atoms, "properties_test.csv", eq_trajObject)
        self.assertTrue(os.path.getsize("properties_test.csv") != 0)      

if __name__ == '__main__':
    tests = [unittest.TestLoader().loadTestsFromTestCase(PropertyCalculationTests)]
    testsuite = unittest.TestSuite(tests)
    result = unittest.TextTestRunner(verbosity=0).run(testsuite)
    sys.exit(not result.wasSuccessful())
