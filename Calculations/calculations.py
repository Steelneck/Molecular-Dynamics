from ase import units
from ase.calculators.eam import EAM
from asap3 import EMT
from ase.data import atomic_masses, atomic_numbers
from ase.md.langevin import Langevin
from ase.lattice.bravais import Lattice
import numpy as np
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.atoms import Atoms
from ase import atoms
from ase.calculators.emt import EMT

# Calculates the specific heat and returns a numpy.float64 with dimensions J/(K*Kg)
# Something wrong with this function. For a bigger crystal the specific heat eneregy is lowered. It should be the opposite.
# Should require more energy to heat up a larger sample. Unclear what causes this. Is it wrong to use the total energy? Is the Langevin function correct? 
# How long should the dyn.run() be going for an accurate results? Reach equilibrium? 
def Specific_Heat(atoms):

    if isinstance(atoms, Lattice) or isinstance(atoms, Atoms):
        # Calculate the bulk mass by 
        mass_per_atom = atomic_masses[atomic_numbers[(atoms.get_chemical_symbols())[0]]]*1.6605402*10**(-27)
        number_of_atoms_in_crystal = len(atoms.numbers)
        bulk_mass = mass_per_atom*number_of_atoms_in_crystal

        # Describe the interatomic interactions with the Effective Medium Theory
        atoms.set_calculator(EMT())

        # Set the momenta corresponding to T=300K
        MaxwellBoltzmannDistribution(atoms, 300 * units.kB)

        # We want to run MD with constant temperture using the Langevin algorithm.
        dyn = Langevin(atoms, 5 * units.fs, units.kB * 300, 0.002)

        temp_vec = np.array([])
        eng_vec = np.array([])

        # Now run the dynamics
        for i in range(2):
            dyn.run(100)
            temp_vec = np.append(temp_vec, (atoms.get_kinetic_energy() / len(atoms)) / (1.5 * units.kB))
            eng_vec = np.append(eng_vec, (atoms.get_potential_energy() / len(atoms)) + (atoms.get_kinetic_energy() / len(atoms)))
        heat_capcity = ((eng_vec[1] - eng_vec[0])*(1.6021765*10**(-19)))/(temp_vec[1] - temp_vec[0]) / bulk_mass
        return heat_capcity
    else:
        return None

