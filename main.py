from ase.lattice.cubic import FaceCenteredCubic
from calculations import Specific_Heat

atoms = FaceCenteredCubic(directions=[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
                          symbol='Al',
                          size=(1, 1, 1),
                          pbc=True)


Specific_Heat(atoms)