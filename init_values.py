""" Initiation of variables and system  """

from ase.lattice.cubic import FaceCenteredCubic
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.md.verlet import VelocityVerlet
from ase import units
from asap3 import EMT

# Set variables
Directions = [[1, 0, 0], [0, 1, 0], [0, 0, 1]] # Orientation of lattice (Miller indices basis vectors of supercell)
Size_X = 10 # How many times fundamental repeat unit is repeated
Size_Y = 10
Size_Z = 10
Symbol = "Cu" # Element specified by atomic symbol e.g. Cu for copper
Pbc = True # Set periodic boundary condition to True or False

# Set up a crystal
atoms = FaceCenteredCubic(directions=Directions,
                          symbol=Symbol,
                          size=(Size_X, Size_Y, Size_Z),
                          pbc=Pbc)

# Set the momenta corresponding to T=300K
MaxwellBoltzmannDistribution(atoms, 300 * units.kB)

# Describe the interatomic interactions with the Effective Medium Theory
atoms.calc = EMT()

