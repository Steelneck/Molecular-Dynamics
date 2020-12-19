import json
import time

from ase import *

"""
Function that will go through all the important input and output from the simulation run in main and saves it too a simulation-specific .json-file.
"""

def write_simulation_to_json(atomName, database, potential, integrator, myAtoms, temperature, MSD, D, L, SHC, Int_T, E_coh, Int_P, Bulk_modulus, Lattice_constant, run_time):
    data = {
        "Simulation input": [
            {
                "Chemical formula": atomName,
                "Database" : database,
                "Potential" : potential,
                "Integrator" : integrator,
                "Lattice structure": str(myAtoms.get_cell().get_bravais_lattice()),
                "a": myAtoms.get_cell_lengths_and_angles()[0],
                "b": myAtoms.get_cell_lengths_and_angles()[1],
                "c": myAtoms.get_cell_lengths_and_angles()[2],
                "Alpha": myAtoms.get_cell_lengths_and_angles()[3],
                "Beta": myAtoms.get_cell_lengths_and_angles()[4],
                "Gamma": myAtoms.get_cell_lengths_and_angles()[5],
                "Volume": myAtoms.get_volume(),
                "Initial Temperature": temperature,
                "Run time": run_time
            }
        ],
        "Simulation output": [
            {
                "Mean Square Displacement": MSD,
                "Self diffusion coefficient": D,
                "Lindemann Criterion": L,
                "Specific Heat": SHC,
                "Internal Temperature": Int_T,
                "Cohesive Energy": E_coh,
                "Internal Pressure": Int_P,
                "Bulk_modulus": Bulk_modulus,
                "Lattice constant": Lattice_constant
            }
        ]
    }
    #Included timestamp (in seconds) in filename to make sure it is always saved as a separate file. .
    filename = "Json/simulation_" + myAtoms.get_chemical_formula() + "_" + str(time.time()) + ".json"

    #Make sure all data is saved to the .json-file
    with open(filename, "w") as file:
        json.dump(data, file, indent = 4)


