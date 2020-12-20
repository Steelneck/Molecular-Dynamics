""" Testing to write atoms object to file """

from ase import *
from ase.lattice.cubic import *
from ase.io import *
from ase.io.cif import read_cif

atoms = FaceCenteredCubic(latticeconstant=None,
                          directions=[[1,0,0],[0,1,0],[0,0,1]],
                          symbol="Cu",
                          size=(10,10,10),
                          pbc=True)

write('testingwrite.cif',
      atoms,
      parallel = True,
      append = False)

cf = read_cif('testingwrite.cif', index=0)

print(next(cf))

print(cf)