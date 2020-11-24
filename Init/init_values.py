""" Initiation of variables and system  """

# 6 of the 7 lattice systems (rhombohedral is not available)
from ase.lattice.cubic import *
from ase.lattice.tetragonal import *
from ase.lattice.orthorhombic import *
from ase.lattice.monoclinic import *
from ase.lattice.triclinic import *
from ase.lattice.hexagonal import *

# Algorithms and calculators for the simulation
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.md.verlet import VelocityVerlet
from ase import units
from asap3 import EMT
from ase.calculators.kim.kim import KIM

# Initiation functions to separate them from variables
from .init_functions import set_lattice
from .init_functions import set_lattice_const

# Pygen 
from pymatgen.ext.matproj import MPRester

m = MPRester('rXy9SNuvaCUyoVmTDjDT')
data = m.query(criteria={"elements": ["Cu"]}, properties=["pretty_formula", "initial_structure"])
# print(((data[0])["spacegroup"])["crystal_system"])
print(type(data))

data1 = (data[0])["initial_structure"]


""" This section will initialize the system based on the variables """
def init_MP():
    i = 0
    atoms_list = []
    while i < len(data): 
        # Set variables for your simulation
        # OBS! Combination of Directions and Miller only works when complete and consistent
        Directions = [[1, 0, 0], [0, 1, 0], [0, 0, 1]] # Orientation of lattice
        Miller = [None, None, None] # Basis of supercell and / or three surfaces
        Size_X = 10 # How many times fundamental repeat unit is repeated
        Size_Y = 10
        Size_Z = 10
        Symbol = str((data[i])["pretty_formula"]) # Element specified by atomic symbol e.g. Cu for copper (OBS! requires string)
        Pbc = (True, True, True) # Set periodic boundary condition to True or False. 
        # if str(((data[i])["spacegroup"])["crystal_system"]) == "cubic":
        #     Bravais = SimpleCubic # Set the lattice
        # elif str(((data[i])["spacegroup"])["crystal_system"]) == "tetragonal":
        #     Bravais = SimpleTetragonal # Set the lattice
        # elif str(((data[i])["spacegroup"])["crystal_system"]) == "monoclinic":
        #     Bravais = SimpleMonoclinic # Set the lattice
        # else:
        #     continue
        Bravais = FaceCenteredCubic
        lc_a = 0 # When lattice constants are zero => FaceCenteredCubic retrieves lc_a from ase
        lc_b = 0
        lc_c = 0
        lc_alpha = 0 # Degrees
        lc_beta = 0
        lc_gamma = 0
        Temperature = 300
        Calculator = EMT()
        i = i + 1

        """ This section is where the user changes values """

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
            For standard Lennard-Jones potential use : LJ_ElliottAkerson_2015_Universal__MO_959249795837_003
        """

        # Set the lattice constant
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


