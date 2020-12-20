"""Demonstrates molecular dynamics with constant energy."""

import ase
from ase.lattice.cubic import FaceCenteredCubic
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.md.verlet import VelocityVerlet
from ase import units
from asap3 import EMT
from asap3 import OpenKIMcalculator


size = 10

# Set up a crystal
#atoms = FaceCenteredCubic(directions=[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
#                          symbol="Cu",
#                          size=(size, size, size),
#                          pbc=True)

atoms = ase.io.read("Fe3Cu_mp-1184260_computed.cif")
atoms = atoms*(size,size,size)
#ase.visualize.view(atoms)

# Describe the interatomic interactions with the Effective Medium Theory
# atoms.calc = EMT()
atoms.set_calculator(OpenKIMcalculator('LJ_ElliottAkerson_2015_Universal__MO_959249795837_003'))

# Set the momenta corresponding to T=300K
MaxwellBoltzmannDistribution(atoms, 300 * units.kB)

# We want to run MD with constant energy using the VelocityVerlet algorithm.
dyn = VelocityVerlet(atoms, 5 * units.fs)  # 5 fs time step.


def printenergy(a=atoms):  # store a reference to atoms in the definition.
    """Function to print the potential, kinetic and total energy."""
    epot = a.get_potential_energy() / len(a)
    ekin = a.get_kinetic_energy() / len(a)
    print('Energy per atom: Epot = %.3feV  Ekin = %.3feV (T=%3.0fK)  '
          'Etot = %.3feV' % (epot, ekin, ekin / (1.5 * units.kB), epot + ekin))

# Now run the dynamics
dyn.attach(printenergy, interval=10)
printenergy()
dyn.run(200)
