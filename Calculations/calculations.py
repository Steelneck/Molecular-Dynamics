from ase import units
from ase.data import atomic_masses, atomic_numbers
from ase.md.langevin import Langevin
import numpy as np
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.calculators.kim.kim import KIM
from asap3 import Trajectory, FullNeighborList

def calc_internal_temperature(myAtoms, trajectoryFileName, timeStepIndex):
    """ Returns the average temperature within parameters """
    
    for i in range(1,timeStepIndex):
        sumTemp = sum(myAtoms.get_temperature())     
        eqTemperature += sumTemp/len(sumTemp)               # Iteration gives the sum of a number of average temperatures

    avgTemperature = eqTemperature/timeStepIndex            # Sampling this gives the average temperature for the atoms
    print("Average internal temperature:", avgTemperature)  
    return(avgTemperature)