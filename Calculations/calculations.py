from ase import units
from ase.data import atomic_masses, atomic_numbers
from ase.md.langevin import Langevin
import numpy as np
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.calculators.kim.kim import KIM
from asap3 import Trajectory, FullNeighborList

"""Function that takes all the atoms-objects after the system reaches equilibrium  (constant total energy, volume and pressure) and writes them over to a new .traj-file. Goes through trajectoryFileName and writes too eq_trajectoryFileName. Uses SuperCellSize to calculate volume."""
def eq_traj(myAtoms, trajectoryFileName, eq_trajectoryFileName, superCellSize):
    traj_non_eq = Trajectory(trajectoryFileName) #Opens the simulation trajectory file
    traj_eq = Trajectory(eq_trajectoryFileName, "w", myAtoms) #new equilibrium trajectory file to write atoms in equilibrium too.
    t = 0
    eq_index = 0
    while t < len(traj_non_eq)-3: #Will check equilibrium conditions for atoms object and stop when equilibrium is reached.
        P_tot_diff = 0
        E_tot_diff = 0
        V_diff = 0
        for i in range(3): #Sums up offset of energy, pressure and volume between timesteps for three different values of i.
            curr_element = t+i
            next_element = t+i+1
            E_tot_1 = (traj_non_eq[curr_element].get_potential_energy() +
                       traj_non_eq[curr_element].get_kinetic_energy()) #Total energy of two consecutive timesteps
            E_tot_2 = (traj_non_eq[next_element].get_potential_energy() +
                       traj_non_eq[next_element].get_kinetic_energy()) 
            P_inst_1 = calc_instantaneous_pressure(myAtoms, trajectoryFileName, curr_element) #instantaneous pressure at two consecutive timesteps
            P_inst_2 = calc_instantaneous_pressure(myAtoms, trajectoryFileName, next_element) 
            V_1 = traj_non_eq[curr_element].get_volume() * superCellSize #Volume of the cell at time two consecutive timesteps
            V_2 = traj_non_eq[next_element].get_volume() * superCellSize 
            V_diff += abs(V_2 - V_1) #Offset in volume between timesteps.
            P_tot_diff += abs(P_inst_2 - P_inst_1) #Offset in pressure
            E_tot_diff += abs(E_tot_2 - E_tot_1) #Offset in total energy
            
        P_tot_diff_mean = P_tot_diff/3 #Mean values of three iterations
        V_diff_mean = V_diff/3
        E_tot_diff_mean = E_tot_diff/3
        
        if E_tot_diff_mean < 0.1 and V_diff_mean < 0.1: #Criteria for equilibrium.
            eq_index = t #saves index of first atom that has reached equilibrium.
            break
        t += 1
    t = len(traj_non_eq) - 1
    while t > eq_index: #while loop that goes through all atoms objects in equilibrium and writes them to new .traj-file
        traj_eq.write(traj_non_eq[t])
        t -= 1
    
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

"""Function to calculate and print the time average of the mean square displacement (MSD) at time t."""
def MSD_calc(myAtoms, trajObject, timeStepIndex):
    try:
        myAtoms.get_masses() # Tries if the attribute exists, skips except if it does
    except AttributeError:
        print("You have not entered a valid system.") # Message for user if no attribute
        return False # Ends the function
    
    time = len(trajObject)-timeStepIndex
    pos_eq = trajObject[-1].get_positions() #position of atoms when system has reached equilibrium
    pos_t = trajObject[1].get_positions() #position of atoms at time t
    diff = pos_t - pos_eq #displacement of all atoms from equilibrium to time t as a vector
    diff_sq = np.absolute(diff)**2 
    MSD = np.sum(diff_sq)/len(myAtoms) #Time averaged mean square displacement.
    print("MSD = ", MSD)
    return(MSD)

def Self_diffuse(trajObject, MSD, timeStepIndex):
    """Function that  calculates the self-diffusion coefficient (D) at time t, based on the value of the mean square displacement."""
    time = len(trajObject)-timeStepIndex #T
    D = MSD/(6*time) #How to connect mean squre displacement to self-diffusion coefficient.
    print("D = ", D)
    return(D)

"""Function that checks the Lindemann criterion which determines if the system is melting or not."""
def Lindemann(trajObject, MSD, timeStepIndex):
     nblist = FullNeighborList(3.5, trajObject[-1]).get_neighbors(1, -1) #Returns 3 lists containing information about nearest neighbors. 3rd list is the square of the distance to the neighbors.
     d = np.sqrt(np.amin(nblist[2])) #distance to the nearest neighbor. Takes the minimum value of nblist.
     L = MSD/d #Lindemann criterion. Expect melting when L>0.1
     if L > 0.1:
         print("MELTING!")
         return True 
     else:
         print("NOT MELTING!")
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
    eqForces = myAtoms.calc.get_forces(atomsAtEquilibrium)    # Use calculator object to get the forces on the atoms. 
    
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
