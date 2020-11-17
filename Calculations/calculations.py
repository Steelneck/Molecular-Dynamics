from ase import units
from ase.data import atomic_masses, atomic_numbers
from ase.md.langevin import Langevin
import numpy as np
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.calculators.kim.kim import KIM

# Calculates the specific heat and returns a numpy.float64 with dimensions J/(K*Kg)
def Specific_Heat(atoms):
    
    try:
        atoms.get_masses() # Tries if the attribute exists, skips except if it does
    except AttributeError:
        print("You have not entered a valid system.") # Message for user if no attribute
        return False # Ends the function

    # Calculates the bulk mass by taking the sum of all the atom masses. 
    bulk_mass=sum(atoms.get_masses())*1.6605402*10**(-27) # converts from atomic mass units to kg

    # Uses the lennard jones potential through openKIM
    calc = KIM("LJ_ElliottAkerson_2015_Universal__MO_959249795837_003")
    atoms.calc = calc

    # Set the momenta corresponding to T=300K
    MaxwellBoltzmannDistribution(atoms, 300 * units.kB)

    # We want to run MD with constant temperture using the Langevin algorithm.
    dyn = Langevin(atoms, 5 * units.fs, units.kB * 300, 0.002)

    temp_vec = np.array([])
    eng_vec = np.array([])

    #MD run with 10 instances in between. Calculates the heat capcity by taking the difference in energy divided by the difference in temperature and 
    #divide everything with the mass of the crystal. 
    for i in range(2):
        dyn.run(10)
        temp_vec = np.append(temp_vec, (atoms.get_kinetic_energy() / len(atoms)) / (1.5 * units.kB))
        eng_vec = np.append(eng_vec, (atoms.get_potential_energy() / len(atoms)) + (atoms.get_kinetic_energy() / len(atoms)))
    
    heat_capcity = ((eng_vec[1] - eng_vec[0])*(1.6021765*10**(-19)))/(temp_vec[1] - temp_vec[0]) / bulk_mass
    return heat_capcity

