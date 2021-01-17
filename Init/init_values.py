""" Initiation of variables and system  """
import math
import shutil
# 6 of the 7 lattice systems (rhombohedral is not available)
import ase.io
import ase
from ase.lattice.cubic import *
from ase.lattice.tetragonal import *
from ase.lattice.orthorhombic import *
from ase.lattice.monoclinic import *
from ase.lattice.triclinic import *
from ase.lattice.hexagonal import *
from ase.atom import *


# Dependencies to run materials project
from pymatgen.ext.matproj import MPRester

# Initiation functions to separate them from variables
from .init_functions import set_lattice_const, set_lattice 

atoms_list = []

# Check if the potential string is empty, then assign standard Lennard Jones potential
def checkKIMpotential(potential):
    if (potential == "") or (potential == " "):
        return "LJ_ElliottAkerson_2015_Universal__MO_959249795837_003"
    else:
        return potential

# Init for ASE
def init(Symbol, Bravais_lattice, Directions, Miller,
            lc_a,lc_b,lc_c,lc_alpha,lc_beta,lc_gamma):


    Lattice_Const = set_lattice_const(lc_a,
                                    lc_b,
                                    lc_c,
                                    lc_alpha,
                                    lc_beta,
                                    lc_gamma)

    # Set up a crystals
    atoms = set_lattice(getattr(ase.lattice.cubic, Bravais_lattice),
                    Lattice_Const,
                    Directions,
                    Miller,
                    1,
                    1,
                    1,
                    Symbol) 

    atoms_list.append(atoms)
    return atoms_list

# Init for Materials project
def init_MP(Criteria_list, API_Key):
    m = MPRester(API_Key) 
    #Loop that takes out each critera for each query
    for criteria in Criteria_list:
        #If there are no elements in data raise an exception and end program
        data = m.query(criteria, properties=['cif', 'spacegroup', 'pretty_formula'])
        if len(data) != 0:
            for i in range(len(data)):
            
                # Takes out the space group and crystal structure query
                crystal_structure = ((data[i])['spacegroup'])['crystal_system']
                pretty_formula = str((data[i])['pretty_formula'])
            
                # Function that skips the element if it not Cubic
                if crystal_structure != 'cubic':
                    continue
                
                
                #Takes out the CIF information and creates a unique file
                cif_Info=(data[i])["cif"]
                f = open(pretty_formula + ".cif", "w+")
                f.write(cif_Info)
                f.close()
                
                #Creates the atom object from the CIF information
                atoms = ase.io.read(pretty_formula + ".cif")
                
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
