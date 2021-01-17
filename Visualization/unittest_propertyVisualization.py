import sys, unittest, os, numpy

from plot_functions import plot_prop_vs_time
from plot_functions import plot_prop_per_simulation
from plot_functions import plot_prop_all_elements
from plot_functions import hist_prop_per_simulation

"""
The following are unittests for the functions in run_plot.py, both for input and output.
"""

class VisualizationTests(unittest.TestCase):

    def test_prop_vs_time_output(self):
        self.assertTrue(len(plot_prop_vs_time("plot_input_test.json", "properties.csv")) != 0)
        
    def test_prop_vs_time_wrong_input(self):
            
        plot1 = plot_prop_vs_time(None, "properties_test.csv")
        plot2 = plot_prop_vs_time("plot_input_test.json", None)

        self.assertIsNone(plot1)
        self.assertIsNone(plot2)

    def test_prop_per_simulation_output(self):
        self.assertTrue(len(plot_prop_per_simulation("plot_input_test.json", "sim_properties_test.csv")) != 0)
        
    def test_prop_per_simulation_wrong_input(self):
        
        plot1 = plot_prop_per_simulation(None, "sim_properties_test.csv")
        plot2 = plot_prop_per_simulation("plot_input_test.json", None)

        self.assertIsNone(plot1)
        self.assertIsNone(plot2)

    def test_all_elements_output(self):
        self.assertTrue(len(plot_prop_all_elements("plot_input_test.json", "sim_properties_test.csv")) != 0)
        
    def test_all_elements_wrong_input(self):
        
        plot1 = plot_prop_all_elements(None, "sim_properties_test.csv")
        plot2 = plot_prop_all_elements("plot_input_test.json", None)
        
        self.assertIsNone(plot1)
        self.assertIsNone(plot2)

    def test_hist_plot_output(self):
        self.assertTrue(len(hist_prop_per_simulation("plot_input_test.json", "sim_properties_test.csv")) != 0)

    def test_hist_plot_wrong_input(self):

        plot1 = hist_prop_per_simulation(None, "sim_properties_test.csv")
        plot2 = hist_prop_per_simulation("plot_input_test.json", None)

        self.assertIsNone(plot1)
        self.assertIsNone(plot2)
        
if __name__ == '__main__':
    tests = [unittest.TestLoader().loadTestsFromTestCase(VisualizationTests)]
    testsuite = unittest.TestSuite(tests)
    result = unittest.TextTestRunner(verbosity=0).run(testsuite)
    sys.exit(not result.wasSuccessful())

    
