import sys, unittest, os, numpy

from ase.lattice.cubic import FaceCenteredCubic
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.md.verlet import VelocityVerlet
from ase import units

from calculations import Heat_Capcity_NVE
#from calculations import Heat_Capcity_NVT
from calculations import internal_temperature
from calculations import cohesive_energy
from calculations import debye_temperature
from calculations import calc_instantaneous_pressure
from calculations import calc_internal_pressure
from calculations import eq_test
from calculations import MSD_calc
from calculations import Self_diffuse
from calculations import Lindemann
from calculations import calc_lattice_constant_fcc_cubic
from calculations import write_time_evolution_to_csv
from calculations import calc_bulk_modulus

from asap3 import Trajectory, EMT

atoms = FaceCenteredCubic(directions=[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
                                  symbol="Cu",
                                  size=(3, 3, 3),
                                  pbc=True)

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
    
    def test_Heat_Capcity_NVE(self):
        HCNVE1 = Heat_Capcity_NVE(atoms, trajObject, 1)
        HCNVE2 = Heat_Capcity_NVE(1, trajObject, 1)
        HCNVE3 = Heat_Capcity_NVE(atoms, 1, 1)
        HCNVE4 = Heat_Capcity_NVE(atoms, trajObject, None)
        
        self.assertIsInstance(HCNVE1, numpy.float64)
        self.assertIsNone(HCNVE2)
        self.assertIsNone(HCNVE3)
        self.assertIsNone(HCNVE4)


    # def test_Heat_Capcity_NVT(self):
    #     HCNVT1 = Heat_Capcity_NVT(atoms, trajObject, 1)
    #     HCNVT2 = Heat_Capcity_NVT(1, trajObject, 1)
    #     HCNVT3 = Heat_Capcity_NVE(atoms, 1, 1)
    #     HCNVT4 = Heat_Capcity_NVE(atoms, trajObject, None)

    #     self.assertIsInstance(HCNVT1, numpy.float64)
    #     self.assertIsNone(HCNVT2)
    #     self.assertIsNone(HCNVT3)
    #     self.assertIsNone(HCNVT4)

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
        self.assertIsInstance(calc_internal_pressure(atoms, trajObject, 1, 1000), float)

    def test_internal_pressure_wrong_input_argument(self):
        internalPressure1 = calc_internal_pressure(None, trajObject, 1, 1000)
        internalPressure2 = calc_internal_pressure(atoms, None, 1, 1000)
        internalPressure3 = calc_internal_pressure(atoms, trajObject, 1, None)
        
        # Then all should return None
        self.assertIsNone(internalPressure1)
        self.assertIsNone(internalPressure2)
        self.assertIsNone(internalPressure3)

    """Unittest for eq_calc"""

    #eq_traj doesnt use atoms yet, so no point in testing that input
    def test_eq_calc_wrong_input_argument(self):
        eq_test2 = eq_test(atoms, None)

        #All should return None
        self.assertIsNone(eq_test2)

    """Unittests for calculation of mean square displacement"""

    def test_MSD_calc_return_type(self):
      	self.assertIsInstance(MSD_calc(atoms, trajObject, -1, 1), float)

    #MSD doesnt use the time input yet so no point in testing it    
    def test_MSD_calc_wrong_input_argument(self):   
        MSD1 = MSD_calc(None, trajObject, -1, 1)
        MSD2 = MSD_calc(atoms, None, -1, 1)
        MSD3 = MSD_calc(atoms, trajObject, None, 1)
        MSD4 = MSD_calc(atoms, trajObject, -1, None)

        #All should return None
        self.assertIsNone(MSD1)
        self.assertIsNone(MSD2)
        self.assertIsNone(MSD3)
        self.assertIsNone(MSD4)

    """Unittests for calculation of Self diffusion coefficient"""

    def test_self_diffuse_return_type(self):
        self.assertIsInstance(Self_diffuse(MSD_calc(atoms, trajObject, -1, 1), len(trajObject), 10, 5), float)

    #Self_diffuse doesnt use the time input yet so no point in testing it 
    def test_Self_diffuse_wrong_input_argument(self):
        D1 = Self_diffuse(None, len(trajObject), 10, 5)
        D2 = Self_diffuse(MSD_calc(atoms, trajObject, -1, 1), None, 10, 5)
        D3 = Self_diffuse(MSD_calc(atoms, trajObject, -1, 1), len(trajObject), None, 5)
        D4 = Self_diffuse(MSD_calc(atoms, trajObject, -1, 1), len(trajObject), 10, None)

        #All should return None
        self.assertIsNone(D1)
        self.assertIsNone(D2)
        self.assertIsNone(D3)
        self.assertIsNone(D4)
        
    def test_Lindemann_return_type(self):
        self.assertIsInstance(Lindemann(trajObject, MSD_calc(atoms, trajObject, -1, 1)), float)

    """Unit tests for internal_temperature"""
    # Test for correct data type (float) returned
    def test_internal_temperature(self):
        self.assertIsInstance(internal_temperature(atoms, trajObject, 1), float)

    # Test for non-negative temperature value
    def test_internal_temperature_not_negative(self):
        self.assertGreaterEqual(internal_temperature(atoms, trajObject, 1), 0)

    # Test for wrong input, expected return is None
    def test_internal_temperature_wrong_input_argument(self):
        self.assertIsNone(internal_temperature(None, trajObject, 1))
        self.assertIsNone(internal_temperature(atoms, None, 1))

    # Test for returned value within 25% of ASE calculation
    def test_internal_temperature_reasonable(self):
        ASE = atoms.get_temperature()
        self.assertAlmostEqual(internal_temperature(atoms, trajObject, 1), ASE, delta=ASE/4)

    """Unit tests for cohesive_energy"""
    # Test for correct data type (float) returned
    def test_cohesive_energy(self):
        self.assertIsInstance(cohesive_energy(atoms, trajObject, 1), float)

    # Test for non-negative, non-zero energy value
    def test_cohesive_energy_positive(self):
        self.assertGreater(cohesive_energy(atoms, trajObject, 1), 0)

    # Test for wrong input, expected return is None
    def test_cohesive_energy_wrong_input_argument(self):
        self.assertIsNone(cohesive_energy(None, trajObject, 1))
        self.assertIsNone(cohesive_energy(atoms, None, 1))

    # Test for returned value within 25% of ASE calculation
    def test_cohesive_energy_reasonable(self):
        ASE = atoms.get_potential_energy()/len(atoms)
        self.assertAlmostEqual(cohesive_energy(atoms, trajObject, 1), ASE, delta=ASE/4)

    """Unit tests for debye_temperature"""
    # Test for correct data type (float) returned
    def test_debye_temperature(self):
        self.assertIsInstance(debye_temperature(trajObject, MSD_calc(atoms, trajObject, -1, 1)), float)

    # Test for non-negative temperature value
    def test_debye_temperature_not_negative(self):
        self.assertGreaterEqual(debye_temperature(trajObject, MSD_calc(atoms, trajObject, -1, 1)), 0)

    # Test for wrong input, expected return is None
    def test_debye_temperature_wrong_input_argument(self):
        self.assertIsNone(debye_temperature(None, trajObject))
        self.assertIsNone(debye_temperature(atoms, None))

    #Lindemann doesnt use the time input yet so no point in testing it 
    def test_Lindemann_wrong_input_argument(self):
        L1 =Lindemann(None, MSD_calc(atoms, trajObject, -1, 1))
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

    """Unittests for calculation of Bulk Modulus"""
    def test_bulk_modulus_wrong_input_argument(self):
        e01, v01, B1 = calc_bulk_modulus(None)
        self.assertIsNone(e01)
        self.assertIsNone(v01)
        self.assertIsNone(B1)

    def test_bulk_modulus_output_type(self):
        e0, v0, B = calc_bulk_modulus(atoms)
        self.assertIsInstance(e0, float)
        self.assertIsInstance(v0, float)
        self.assertIsInstance(B, float)
        
    # def test_minimum_energy_bulk_modulus(self):
    #     # Warning increasing lattice constant even more, to 20 * cell will make the function execute but values will be crazy. No error handling for this yet.
    #     cell = atoms.get_cell()
    #     atoms.set_cell(cell * 10, scale_atoms=True)     # Modify cell an absurd amount. EOS will give error here
    #     e0, v0, B = calc_bulk_modulus(atoms)
    #     self.assertIsNone(e0)
    #     self.assertIsNone(v0)
    #     self.assertIsNone(B)
    #     atoms.set_cell(cell, scale_atoms=True)          # Reset cell


        # atoms.set_cell(cell * 6, scale_atoms=True)     # Modify cell an absurd amount. Minimum at ends gives error here.
        # e0, v0, B = calc_bulk_modulus(atoms)
        # self.assertIsNone(e0)
        # self.assertIsNone(v0)
        # self.assertIsNone(B)
        # atoms.set_cell(cell, scale_atoms=True)          # Reset cell

    
    def test_csv_writer_wrong_input_argument(self):
        csv1 = write_time_evolution_to_csv(None, "properties_test.csv", trajObject, 1, 10, 5)
        csv2 = write_time_evolution_to_csv(atoms, None, trajObject, 1, 10, 5)
        csv3 = write_time_evolution_to_csv(atoms, "properties_test.csv", None, 1, 10, 5)
        csv4 = write_time_evolution_to_csv(atoms, "properties_test.csv", trajObject, None, 10, 5)
        csv5 = write_time_evolution_to_csv(atoms, "properties_test.csv", trajObject, 1, None, 5)

        self.assertIsNone(csv1)
        self.assertIsNone(csv2)
        self.assertIsNone(csv3)
        self.assertIsNone(csv4)
        self.assertIsNone(csv5)
        
        
    def test_csv_writer_check_csv(self):
    #Check that the .csv-file exists and is not empty.
         write_time_evolution_to_csv(atoms, "properties_test.csv", trajObject, 1, 10, 5)
         self.assertTrue(os.path.getsize("properties_test.csv") != 0)      

if __name__ == '__main__':
    tests = [unittest.TestLoader().loadTestsFromTestCase(PropertyCalculationTests)]
    testsuite = unittest.TestSuite(tests)
    result = unittest.TextTestRunner(verbosity=0).run(testsuite)
    sys.exit(not result.wasSuccessful())

