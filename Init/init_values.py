""" Initiation of variables and system  """
import math
import shutil


from ase.atom import *
import ase.io

#from asap3 import OpenKIMcalculator 
from asap3 import Trajectory

# Algorithms and calculators for the simulation
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution, Stationary
from ase.md.verlet import VelocityVerlet
import asap3.md.verlet
from ase.md.langevin import Langevin
from ase import units
from asap3 import EMT
import asap3
from ase.calculators.kim.kim import KIM

# Initiation functions to separate them from variables
from .init_functions import create_vacancy, find_crystal_center, set_lattice
from .init_functions import set_lattice_const

# Dependencies to run materials project
from pymatgen.ext.matproj import MPRester
from pymatgen.io.cif import CifParser
from .init_functions import insert_impurity
import Init.init_functions

input_list = []

# Check if the potential string is empty, then assign standard Lennard Jones potential
# def checkKIMpotential(potential):
#     if (potential == "") or (potential == " "):
#         return "LJ_ElliottAkerson_2015_Universal__MO_959249795837_003"
#     else:
#         return potential
    


# Init for Materials project
def init_MP(EMT_Check,openKIM_Check,Verlocity_Verlet_Check,KIM_potential,Critera_list,
                                Vacancy, Impurity, Impurity_ele, Temperature,
                                Size_X,Size_Y,Size_Z,API_Key,PBC):
    #API key to fetch data from Materials project
    m = MPRester(API_Key) 
    
    #Loop that takes out each critera for each query
    for criteria in Critera_list:

        #If there are no elements in data raise an exception and end program
        data = m.query(criteria, properties=['cif', 'spacegroup', 'pretty_formula'])
        if len(data) != 0:
            for i in range(len(data)):
            
                # Takes out the space group and crystal structure query
                space_group = ((data[i])['spacegroup'])['symbol']
                crystal_structure = ((data[i])['spacegroup'])['crystal_system']
                pretty_formula = str((data[i])['pretty_formula'])
            
                # Function that skips the element if it not an FCC crystal
                if space_group[0] != 'F' or crystal_structure != 'cubic':
                    continue
                
                #Takes out the CIF information and creates a unique file
                cif_Info=(data[i])["cif"]
                f = open(pretty_formula + ".cif", "w+")
                f.write(cif_Info)
                f.close()
                
                #Creates the atom object from the CIF information
                atoms = ase.io.read(pretty_formula + ".cif")
                atoms.set_pbc(PBC)
                
                #Creates a supercell
                atoms = atoms*(Size_X,Size_Y,Size_Z)

                #places impurity in the crystal 
                if Impurity == True:
                    atom_pos = find_crystal_center(atoms) # Returns a center position in the crystal
                    insert_impurity(atoms, Impurity_ele, atom_pos) # Insert "foregin" atom in the crystal

                #Places vacancy in the crytal
                if Vacancy == True:
                    create_vacancy(atoms) # Create a vacancy

                # Set the momenta corresponding to desired temperature when running Verlocity Verlet
                if Verlocity_Verlet_Check == True:
                    MaxwellBoltzmannDistribution(atoms, Temperature * units.kB)
                    Stationary(atoms) # Set linear momentum to zero

                # Interatomic potential
                if (EMT_Check == True) and (openKIM_Check == False):
                    atoms.calc = EMT()
                elif (EMT_Check == False) and (openKIM_Check == True):
                    #Sets the potential for openKIM. If none is given returns standard Lennard-Jones
                    potential = checkKIMpotential(KIM_potential)
                    atoms.calc = KIM(potential, options={"ase_neigh": True})
                    #atoms.set_calculator(OpenKIMcalculator(potential))
                else:
                    raise Exception("EMT=openKIM. Both cannot be true/false at the same time!")
                
                #Moves the trajectory file to another folder after it has been used
                shutil.move(pretty_formula + ".cif", "CIF/" + pretty_formula + ".cif")
                atoms_list.append(atoms)
        else:
            raise Exception("The query/queries returned no elements!") 
    
    # If atoms_list is empty raise an expetion and end program
    if len(atoms_list) == 0:
        raise Exception("The query/queries contains no FCC crystals")
    else:
        return atoms_list
