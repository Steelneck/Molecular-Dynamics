""" Initiation of variables and system  """

import math

# 6 of the 7 lattice systems (rhombohedral is not available)
from ase.lattice.cubic import *
from ase.lattice.tetragonal import *
from ase.lattice.orthorhombic import *
from ase.lattice.monoclinic import *
from ase.lattice.triclinic import *
from ase.lattice.hexagonal import *
from ase import Atoms

# Algorithms and calculators for the simulation
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.md.verlet import VelocityVerlet
from ase import units
from asap3 import EMT
from ase.calculators.kim.kim import KIM

# Initiation functions to separate them from variables
from .init_functions import set_lattice
from .init_functions import set_lattice_const
from .init_functions import from_dictionary_to_atoms

# Dependencies to run materials project
from pymatgen.ext.matproj import MPRester
from pymatgen.io.cif import CifParser

""" This section is where the user changes values """

# Set variables for your simulation
# OBS! Combination of Directions and Miller only works when complete and consistent
Directions = [[1, 0, 0], [0, 1, 0], [0, 0, 1]] # Orientation of lattice
Miller = [None, None, None] # Basis of supercell and / or three surfaces
Size_X = 10 # How many times fundamental repeat unit is repeated
Size_Y = 10
Size_Z = 10
Symbol = "Cu" # Element specified by atomic symbol e.g. Cu for copper (OBS! requires string)
Pbc = (True, True, True) # Set periodic boundary condition to True or False. 
Bravais = FaceCenteredCubic # Set the lattice
lc_a = 0 # When lattice constants are zero => FaceCenteredCubic retrieves lc_a from ase
lc_b = 0
lc_c = 0
lc_alpha = 0 # Degrees
lc_beta = 0
lc_gamma = 0
Temperature = 300
Calculator = EMT()
steps = 1000 # Timesteps for dyn.run
interval = 10 # Writes in traj at n timestep


""" The following Bravais lattices can be used:
 SimpleCubic                 Lattice constant: a
 FaceCenteredCubic           Lattice constant: a (set a=0 to let ase set the constant)
 BodyCenteredCubic           Lattice constant: a
 SimpleTetragonal            Lattice constant: a,c
 CenteredTetragonal          Lattice constant: a,c
 SimpleOrthorhombic          Lattice constant: a,b,c
 BaseCenteredOrthorhombic    Lattice constant: a,b,c
 FaceCenteredOrthorhombic    Lattice constant: a,b,c
 BodyCenteredOrthorhombic    Lattice constant: a,b,c
 SimpleMonoclinic            Lattice constant: a,b,c,alpha
 BaseCenteredMonoclinic      Lattice constant: a,b,c,alpha
"""

""" The following Bravais lattices cannot be used:
 Diamond # (requires a basis - check ase how to define new lattice with Factory)
 Triclinic # Not sure how this one works, needs all constants a,b,c,alpha,beta,gamma
 Hexagonal # OBS! Currently broken, see ase/lattice/hexagonal.py line 60
 HexagonalClosedPacked # OBS! Currently broken (and requires a basis)
 Hexagonal / HexagonalClosedPacked needs constant a,c
 Graphite (requires a basis)
"""

""" The following Calculators can be used:
EMT()
    ASAP3 built in effective medium theory. Works for Ni, Cu, Pd, Ag, Pt and Au (and their alloys).
KIM('Insert_openKIM_potential_here')
    openKIM potentials can be found from https://openkim.org/
    For standard lenard-jones potential use: LJ_ElliottAkerson_2015_Universal__MO_959249795837_003
"""

""" This section is only for parameters is concerned with Materials project 
    The size parameters, calculator and temperature above is also used when running materials project
"""

m = MPRester('rXy9SNuvaCUyoVmTDjDT') # Insert your API-Key from https://materialsproject.org/
criteria= {"elements": ["Cu"]}

""" MongoDB query to get desired data from materialsproject 
        Quries to use can be found on https://docs.mongodb.com/manual/reference/operator/query/
Note: Right now this function sorts out everything except for FCC crystals
"""

""" Decide timestepindex for the traj file """

def timestepindex(timesteps, traj_interval):
    trajindexes = math.floor(timesteps/traj_interval)
    return trajindexes

timeStepIndex = timestepindex(steps, interval)
atoms_list = []

def init():
    Lattice_Const = set_lattice_const(lc_a,
                                    lc_b,
                                    lc_c,
                                    lc_alpha,
                                    lc_beta,
                                    lc_gamma)

    # Set up a crystals
    atoms = set_lattice(Bravais,
                    Lattice_Const,
                    Directions,
                    Miller,
                    Size_X,
                    Size_Y,
                    Size_Z,
                    Symbol,
                    Pbc) 

        
    # Set the momenta corresponding to T=300K 
    # (Note: Create a higher order function)
    MaxwellBoltzmannDistribution(atoms, Temperature * units.kB)

    # Describe the interatomic interactions with the Effective Medium Theory
    # (Note: Create a higher ordet function)
    atoms.calc = Calculator

    atoms_list.append(atoms)
    return atoms_list

def init_MP():    
    #If there are no elements in data raise an exception and end program
    data = m.query(criteria, properties=['cif', 'spacegroup', 'pretty_formula'])
    if len(data) != 0:
        for i in range(len(data)):
            
            # Takes out the space group and crystal structure query
            space_group = ((data[i])['spacegroup'])['symbol']
            crystal_structure = ((data[i])['spacegroup'])['crystal_system']
            
            # Function that skips the element if it not an FCC crystal
            if space_group[0] != 'F' or crystal_structure != 'cubic':
                continue
            
            # Ordered dictionary of the CIF
            cif_Info=(CifParser.from_string((data[i])["cif"])).as_dict()

            #Pretty formula for the element
            pretty_formula = str((data[i])['pretty_formula'])

            # Function that returns an atomobject depending on the information from the CIF
            atoms = from_dictionary_to_atoms(cif_Info, pretty_formula, Size_X, Size_Y, Size_Z)

            # Set the momenta corresponding to T=300K 
            # (Note: Create a higher order function)
            MaxwellBoltzmannDistribution(atoms, Temperature * units.kB)

            # Describe the interatomic interactions with the Effective Medium Theory
            # (Note: Create a higher ordet function)
            atoms.calc = Calculator

            atoms_list.append(atoms)

        # If atoms_list is empty raise an expetion and end program
        if len(atoms_list) == 0:
            raise Exception("The query contains no FCC crystals")
        else:
            return atoms_list
    else:
        raise Exception("The query returned no elements!") 