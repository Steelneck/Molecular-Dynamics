import sys, unittest
from ase.lattice.cubic import FaceCenteredCubic
from calculations import Specific_Heat
from calculations import calc_instantaneous_pressure, calc_internal_pressure
from calculations import internal_temperature
import numpy
from asap3 import Trajectory, EMT

atoms = FaceCenteredCubic(directions=[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
                                  symbol="Cu",
                                  size=(3, 3, 3),
                                  pbc=True)

not_bravice_lattice = 1

# Create trajectory object
trajWriter = Trajectory("test.traj", "w", atoms)
trajWriter.write(atoms)
trajWriter.write(atoms)
trajWriter.write(atoms)
trajWriter.write(atoms)
trajObject = Trajectory("test.traj")

# Attach calculator to atoms object
atoms.calc = EMT()

class PropertyCalculationTests(unittest.TestCase):
 

    def test_specific_heat(self):
        self.assertIsInstance(Specific_Heat(atoms), numpy.float64)
        
    def test_specific_heat_not_bravice_lattice(self):
        self.assertFalse(Specific_Heat(not_bravice_lattice))

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

    def test_internal_temperature(self):
        self.assertIsInstance(internal_temperature(atoms, trajObject, 1000), float)

    def test_internal_temperature_not_negative(self):
        self.assertGreaterEqual(internal_temperature(atoms, trajObject, 1000), 0)


if __name__ == '__main__':
    tests = [unittest.TestLoader().loadTestsFromTestCase(PropertyCalculationTests)]
    testsuite = unittest.TestSuite(tests)
    result = unittest.TextTestRunner(verbosity=0).run(testsuite)
    sys.exit(not result.wasSuccessful())