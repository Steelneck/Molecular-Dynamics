import sys, unittest, os, numpy
from main import simulation
import json

"""
simulation(EMT_Check,openKIM_Check, Lennard_Jones_Check, LJ_epsilon,
                        LJ_sigma, LJ_cutoff, Verlocity_Verlet_Check, 
                        Langevin_Check, Langevin_friction, time_step, KIM_potential,
                        ASE, Symbol, Materials_project,API_Key,Criteria_list, 
                        Vacancy, Impurity, Impurity_ele,
                        Temperature, Steps, Interval,Size_X, Size_Y, Size_Z,
                        PBC, Bravais_lattice,Directions,Miller,
                        lc_a,lc_b,lc_c,lc_alpha,lc_beta,lc_gamma,run_Optimade,Optimade_name,Optimzed_volume)
"""

""" 
    To run this test one need to remove all the stored .json files in the folder Json
    OBS: Dont forget to store your .json files!
    Do this command (if you are in the main folder) before doing the test: rm Json/*.json
"""

""" Tests if the system raises expetion with feeding it wrong input """
class acceptanceTests(unittest.TestCase):

    def test_system_is_fed_wrong(self):
        with self.assertRaises(Exception): 
            simulation(True,True, False, [0.010323], [3.40], 6.625, True, 
                         False, 0.002, 5, "", True, "Cu", False,"rXy9SNuvaCUyoVmTDjDT",[{"elements" : ["Cu"]}], 
                        False, False, "Ni", 300, 10000, 100,5, 5, 5,
                        [True,True,True], "FaceCenteredCubic", [[1, 0, 0], [0, 1, 0],
                        [0, 0, 1]], [None, None, None],5.256,0,0,0,0,0,False,"run1", False)
        with self.assertRaises(Exception): 
            simulation(True,False, True, [0.010323], [3.40], 6.625, True, 
                        False, 0.002, 5, "", True, "Cu", False,"rXy9SNuvaCUyoVmTDjDT",[{"elements" : ["Cu"]}], 
                        False, False, "Ni", 300, 10000, 100,5, 5, 5,
                        [True,True,True], "FaceCenteredCubic", [[1, 0, 0], [0, 1, 0],
                        [0, 0, 1]], [None, None, None],5.256,0,0,0,0,0,False,"run1", False)
        with self.assertRaises(Exception): 
            simulation(True,False, False, [0.010323], [3.40], 6.625, True, 
                        True, 0.002, 5, "", True, "Cu", False,"rXy9SNuvaCUyoVmTDjDT",[{"elements" : ["Cu"]}], 
                        False, False, "Ni", 300, 10000, 100,5, 5, 5,
                        [True,True,True], "FaceCenteredCubic", [[1, 0, 0], [0, 1, 0],
                        [0, 0, 1]], [None, None, None],5.256,0,0,0,0,0,False,"run1", False)
        with self.assertRaises(Exception): 
            simulation(True,False, False, [0.010323], [3.40], 6.625, True, 
                        False, 0.002, 5, "", True, "Cu", True,"rXy9SNuvaCUyoVmTDjDT",[{"elements" : ["Cu"]}], 
                        False, False, "Ni", 300, 10000, 100,5, 5, 5,
                        [True,True,True], "FaceCenteredCubic", [[1, 0, 0], [0, 1, 0],
                        [0, 0, 1]], [None, None, None],5.256,0,0,0,0,0,False,"run1", False)

    """ Sees that an FCC run with copper calculated with EMT() returns float numbers """  
    def test_system_returns_values(self):
        simulation(True,False, False, [0.010323], [3.40], 6.625, True, 
                        False, 0.002, 5, "", True, "Cu", False,"rXy9SNuvaCUyoVmTDjDT",[{"elements" : ["Cu"]}], 
                        False, False, "Ni", 300, 1000, 10, 5, 5, 5,
                        [True,True,True], "FaceCenteredCubic", [[1, 0, 0], [0, 1, 0],
                        [0, 0, 1]], [None, None, None],5.256,0,0,0,0,0,False,"run1", False)


        filename = "Json/simulation_Cu"

        for f_name in os.listdir('Json'):
            if f_name.startswith(filename):
                with open(f_name) as json_file:
                    values = json.load(json_file)                        

                    self.assertIsInstance((values["Simulation output"])["Mean Square Displacement"], numpy.float64)
                    self.assertIsInstance((values["Simulation output"])["Self diffusion coefficient"], numpy.float64)
                    self.assertIsInstance((values["Simulation output"])["Lindemann Criterion"], numpy.float64)
                    self.assertIsInstance((values["Simulation output"])["Specific Heat"], numpy.float64)
                    self.assertIsInstance((values["Simulation output"])["Internal Temperature"], numpy.float64)
                    self.assertIsInstance((values["Simulation output"])["Cohesive Energy"], numpy.float64)
                    self.assertIsInstance((values["Simulation output"])["Internal Pressure"], numpy.float64)
                    self.assertIsInstance((values["Simulation output"])["Bulk_modulus"], numpy.float64)
                    self.assertIsInstance((values["Simulation output"])["Lattice constant"], numpy.float64)

    """ Tests that values for Argon is accurate"""
    def test_system_argon_values(self):
        
        simulation(False,False, True, [0.010323],
                [3.40], 6.625, True, False, 0.002, 1, "",
                True, "Ar", False ,"rXy9SNuvaCUyoVmTDjDT", [{"elements" : ["Cu"]}], 
                False, False, "Cu", 300, 20000, 100, 15, 15, 15,
                [True,True,True], "FaceCenteredCubic", [[1, 0, 0], [0, 1, 0],
                [0, 0, 1]], [None, None, None],5.256,0,0,0,0,0,False,"run1", False)
    
        filename = "Json/simulation_Ar"

        for f_name in os.listdir('Json'):
            if f_name.startswith(filename):
                with open(f_name) as json_file:
                    values = json.load(json_file)


                    self.assertAlmostEqual((values["Simulation output"])["Specific Heat"], 312)
                    self.assertAlmostEqual((values["Simulation output"])[["Cohesive Energy"]], 0.08)
                    self.assertAlmostEqual((values["Simulation output"])[["Bulk_modulus"]], 6.1)
                    self.assertAlmostEqual((values["Simulation output"])[["Internal Pressure"]], 0)


if __name__ == '__main__':
    tests = [unittest.TestLoader().loadTestsFromTestCase(acceptanceTests)]
    testsuite = unittest.TestSuite(tests)
    result = unittest.TextTestRunner(verbosity=0).run(testsuite)
    sys.exit(not result.wasSuccessful())