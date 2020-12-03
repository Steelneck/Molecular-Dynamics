""" Initiation of variables and system  """
from ase import Atoms
from ase.spacegroup import crystal
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
# Ni, Cu, Pd, Ag, Pt, and Au,
data = m.query(criteria={"elements": ["Cu"]}, properties=["cif", 'spacegroup', 'pretty_formula'])
print(len(data))
Temperature = 300
Calculator = EMT()
size_x = 10
size_y = 10
size_z = 10

""" This section will initialize the system based on the variables """
def init_MP():
    atoms_list = []
    for i in range(len(data)):
        
        pretty_formula = str((data[i])['pretty_formula'])
        space_group = ((data[i])['spacegroup'])['symbol']
        crystal_structure = ((data[i])['spacegroup'])['crystal_system']
        
        
        if space_group[0] != 'F' or crystal_structure != 'cubic':
           print('Not FCC!')
           continue

        cif_Info=(CifParser.from_string((data[i])["cif"])).as_dict()
        print(cif_Info)

        atoms = from_dictionary_to_atoms(cif_Info, pretty_formula)

        # Set the momenta corresponding to T=300K 
        # (Note: Create a higher order function)
        MaxwellBoltzmannDistribution(atoms, Temperature * units.kB)

        # Describe the interatomic interactions with the Effective Medium Theory
        # (Note: Create a higher ordet function)
        atoms.calc = Calculator

        atoms_list.append(atoms)
    return atoms_list

def from_dictionary_to_atoms(dictionary, symbol):

    # Use space group instead of atoms? Would that be better? 
    chemical_formula = str((dictionary[symbol])['_chemical_formula_sum'])
    a = float((dictionary[symbol])['_cell_length_a'])
    b = float((dictionary[symbol])['_cell_length_b'])
    c = float((dictionary[symbol])['_cell_length_c'])
    alpha = float((dictionary[symbol])['_cell_angle_alpha'])
    beta = float((dictionary[symbol])['_cell_angle_beta'])
    gamma = float((dictionary[symbol])['_cell_angle_gamma'])
    nr_of_atoms = len(((dictionary[symbol])['_atom_site_fract_x']))

    position = []
    for i in range(nr_of_atoms):
        post = list()
        x = float(((dictionary[symbol])['_atom_site_fract_x'])[i])*a
        y = float(((dictionary[symbol])['_atom_site_fract_y'])[i])*b
        z = float(((dictionary[symbol])['_atom_site_fract_z'])[i])*c
        pos = (x,y,z)
        position.append(pos)
    print(position)        

    atoms =Atoms(symbols= chemical_formula,
                positions=position,
                cell=[a, b, c, alpha, beta, gamma])

    atoms = atoms*(size_x,size_y,size_z)
    #a = np.matrix('5 5 5; 5 5 5; 5 5 5')
    #print(a)
    #atoms = make_supercell(atom, a.any())

    return atoms

