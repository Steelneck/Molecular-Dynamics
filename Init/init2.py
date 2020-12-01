""" Initiation of variables and system  """
from ase import Atoms
import numpy as np

# Algorithms and calculators for the simulation
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.md.verlet import VelocityVerlet
from ase import units
from asap3 import EMT
from ase.calculators.kim.kim import KIM
from ase.build import make_supercell

# Pymagen 
from pymatgen.ext.matproj import MPRester
from pymatgen.io.cif import CifParser

# Insert API-Key to be able to run. 
m = MPRester('rXy9SNuvaCUyoVmTDjDT')

#MongoDB query to get data from materials project. Properties must allways show "cif" to work.
data = m.query(criteria={"elements": ["Cu"]}, properties=["cif"])
#print(data)
Temperature = 300
Calculator = EMT()

""" This section will initialize the system based on the variables """
def init_MP():
    i = 0
    atoms_list = []
    while i < 2: 
        cif_Info=(CifParser.from_string((data[i])["cif"])).as_dict()
        
        atoms = from_dictionary_to_atoms(cif_Info)

        #print(cif_Info['Cu'])

        # Set the momenta corresponding to T=300K 
        # (Note: Create a higher order function)
        MaxwellBoltzmannDistribution(atoms, Temperature * units.kB)

        # Describe the interatomic interactions with the Effective Medium Theory
        # (Note: Create a higher ordet function)
        atoms.calc = Calculator

        atoms_list.append(atoms)
        i = i+1
    return atoms_list

def from_dictionary_to_atoms(dictionary):

    # Use space group instead of atoms? Would that be better? 
    symbol = str((dictionary['Cu'])['_chemical_formula_sum'])
    a = float((dictionary['Cu'])['_cell_length_a'])
    b = float((dictionary['Cu'])['_cell_length_b'])
    c = float((dictionary['Cu'])['_cell_length_c'])
    alpha = float((dictionary['Cu'])['_cell_angle_alpha'])
    beta = float((dictionary['Cu'])['_cell_angle_beta'])
    gamma = float((dictionary['Cu'])['_cell_angle_gamma'])
    nr_of_atoms = int((dictionary['Cu'])['_cell_formula_units_Z'])

    #is atom site fract the position? 
    position = []
    for i in range(nr_of_atoms):
        post = list()
        x = float(((dictionary['Cu'])['_atom_site_fract_x'])[i])
        y = float(((dictionary['Cu'])['_atom_site_fract_y'])[i])
        z = float(((dictionary['Cu'])['_atom_site_fract_z'])[i])
        pos = (x,y,z)
        position.append(pos)        

    atoms =Atoms(symbols= symbol,
                positions=position,
                cell=[a, b, c, alpha, beta, gamma])

    atoms = atoms*(10,10,10)
    #a = np.matrix('5 5 5; 5 5 5; 5 5 5')
    #print(a)
    #atoms = make_supercell(atom, a.any())

    return atoms
