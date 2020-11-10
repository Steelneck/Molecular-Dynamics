""" Initiation of variables and system  """

from ase.lattice.cubic import FaceCenteredCubic
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.md.verlet import VelocityVerlet
from ase import units
from asap3 import EMT

Directions = [[1, 0, 0], [0, 1, 0], [0, 0, 1]] # Orientation of lattice (Miller indices basis vectors of supercell)
Size = 10 # How many times fundamental repeat unit is repeated
Symbol = "Cu" # Element specified by atomic symbol e.g. Cu for copper
Pbc = True # Set periodic boundary condition to True or False

# Set up a crystal
atoms = FaceCenteredCubic(directions=Directions,
                          symbol=Symbol,
                          size=(Size, Size, Size),
                          pbc=Pbc)

# Set the momenta corresponding to T=300K
MaxwellBoltzmannDistribution(atoms, 300 * units.kB)

# Describe the interatomic interactions with the Effective Medium Theory
atoms.calc = EMT()
