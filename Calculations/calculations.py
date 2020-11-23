from ase import units
from ase.data import atomic_masses, atomic_numbers
import numpy as np
from asap3 import Trajectory, FullNeighborList
from ase.calculators.eam import EAM
from ase.build import bulk

#Function that takes the atoms-objects that have reached equilibrium and writes them over to a new .traj-file.
def eq_traj(a):
    traj_non_eq = Trajectory("atoms.traj")
    traj_eq = Trajectory("atoms_eq.traj", "w", a) #new equilibrium trajectory file
    n = 1
    while n < len(traj_non_eq): #Goes through all objects in traj_non_eq and checks whether or not they are in equilibrium.
        ediff = abs((traj_non_eq[n].get_potential_energy() / len(a)) -
                    (traj_non_eq[n].get_kinetic_energy() / len(a))) #epot-ekin
        if ediff < 0.003: #Criteria for equilibrium.
            traj_eq.write(traj_non_eq[n]) #write that object to the new .traj-file
        n += 1

# Calculates the specific heat and returns a numpy.float64 with dimensions J/(K*Kg)
def Specific_Heat(a):
    
    try:
        a.get_masses() # Tries if the attribute exists, skips except if it does
    except AttributeError:
        print("You have not entered a valid system.") # Message for user if no attribute
        return False # Ends the function
    traj = Trajectory("atoms_eq.traj")
    bulk_mass=sum(a.get_masses())*1.6605402*10**(-27)
    temp_diff = (traj[1].get_kinetic_energy() /len(traj[1]) - traj[0].get_kinetic_energy() /len(traj[0])) / (1.5 * units.kB)  #Temperature difference between two runs when system has reached equilibrium
    pot_energy_diff = (traj[1].get_potential_energy() /len(traj[1]) 
                        - traj[0].get_potential_energy() /len(traj[0])) # potential energy difference when ystem has reached equilibrium

    kin_energy_diff = (traj[1].get_kinetic_energy() /len(traj[1]) 
                            - traj[0].get_kinetic_energy()/len(traj[0])) # potential energy difference when ystem has reached equilibrium
    
    heat_capcity = abs(((pot_energy_diff + kin_energy_diff)*(1.6021765*10**(-19)))/(temp_diff) / bulk_mass)
    print("C_p = ", heat_capcity, "[J/K*Kg]")
    return heat_capcity

def MSD_calc(a, t):
    try:
        a.get_masses() # Tries if the attribute exists, skips except if it does
    except AttributeError:
        print("You have not entered a valid system.") # Message for user if no attribute
        return False # Ends the function
    """Function to calculate and print the time average of the mean square displacement (MSD) at time t."""
    traj_MSD = Trajectory("atoms_eq.traj")  #time evolution of trajectory
    time = len(traj_MSD)-t
    pos_eq = traj_MSD[-1].get_positions() #position of atoms when system has reached equilibrium
    pos_t = traj_MSD[t].get_positions() #position of atoms at time t
    diff = pos_t - pos_eq #displacement of all atoms from equilibrium to time t as a vector
    diff_sq = np.absolute(diff)**2 
    MSD = np.sum(diff_sq)/len(a) #Time averaged mean square displacement.
    print("MSD = ", MSD)
    return(MSD)

def Self_diffuse(MSD, t):
    """Function that  calculates the self-diffusion coefficient (D) at time t, based on the value of the mean square displacement."""
    traj_MSD = Trajectory("atoms_eq.traj")  
    time = len(traj_MSD)-t
    D = MSD/(6*time) #How to connect mean squre displacement to self-diffusion coefficient.
    print("D = ", D)
    return(D)

def Lindemann(a, MSD):
     nblist = FullNeighborList(3.5, a).get_neighbors(1, -1) #Returns 3 lists containing information about nearest neighbors. 3rd list is the square of the distance to the neighbors.
     d = np.sqrt(np.amin(nblist[2])) #distance to the nearest neighbor. Takes the minimum value of nblist.
     L = MSD/d #Lindemann criterion. Expect melting when L>0.1
     if L > 0.1:
         #print("MELTING!")
         print(L)
         return True
     else:
         #print("Not melting.")
         print(L)
         return False

def calc_instantaneous_pressure(myAtoms, trajectoryFileName, timeStepIndex):
    """ Calculates instantaneous pressure at time timeStepIndex * deltaT, i.e. P(n*deltaT)."""
    traj = Trajectory(trajectoryFileName)
    equilibriumIndex = timeStepIndex
    atomsAtEquilibrium = traj[equilibriumIndex]             # Retrieve atoms some amount of timesteps in when system has reached equilibrium
    eqPos = atomsAtEquilibrium.get_positions()              # Atom 3D coords.
    eqTemperature = atomsAtEquilibrium.get_temperature()    # Temperature at time t = equilibriumIndex * deltaT.
    eqPotEn = atomsAtEquilibrium.get_potential_energy()     # Get potetntial energy from atoms object at time t = equilibriumIndex * deltaT, should be the same for all. 
    N = len(atomsAtEquilibrium)                             # Number of atoms in object
    eqForces = atoms.calc.get_forces(atomsAtEquilibrium)    # Use calculator object to get the forces on the atoms. 
    
    # Instantaneos pressure 
    posForce = np.sum(eqPos * eqForces, axis=1)             # Dot product of all forces and positions
    instantP = (N * units.kB * eqTemperature) / eqPotEn + np.sum(posForce)

    #print(instantP)
    return(instantP)

def calc_internal_pressure(myAtoms, trajectoryFileName, iterations):
    """ Internal pressure is the MD average of the instantaneous pressures """
    M = iterations                                          # M = (iterations * deltaT) / deltaT

    # MD average
    allInstantPressures = 0
    for n in range(1,iterations):
        allInstantPressures += calc_instantaneous_pressure(myAtoms, trajectoryFileName, n)
    
    internalPressure = allInstantPressures / M              # Internal pressure is the MD average of the instantaneous pressures
    print("Internal Pressure", internalPressure)
    return(internalPressure)

def internal_temperature(myAtoms, timeStepIndex):
    """ Returns the average temperature within parameters """
    eqTemp = 0                                             

    for i in range(1, timeStepIndex):                       
        eqTemp += myAtoms.get_temperature()                 # Sum returned value from ASE function over timesteps for sampling
     
    print("Internal temperature:", eqTemp/timeStepIndex, "[K]")  
    return(eqTemp/timeStepIndex)                            # Average over number of samples, return a final value

def cohesive_energy():
    """ Returns the cohesive energy of the system """

    alloy = bulk('Cu', 'fcc', 4.05)

    EAMresult = EAM(potential='Calculations/Cu_Zhou04.eam.alloy')
    EAMresult.write_potential('Calculations/new.eam.alloy')
    alloy.calc = EAMresult
    cohEnergy = alloy.get_potential_energy()

    print("Cohesive energy:", cohEnergy, "[eV]")
    return(cohEnergy)