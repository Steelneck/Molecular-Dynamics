from ase import units
from ase.data import atomic_masses, atomic_numbers
import numpy as np
from asap3 import Trajectory, FullNeighborList
from ase.calculators.eam import EAM
from ase.build import bulk
from ase.io import read

from asap3 import EMT ###!!! Temporary

"""Function that takes all the atoms-objects after the system reaches equilibrium  (constant total energy, volume and pressure) and writes them over to a new .traj-file. Goes through trajectoryFileName and writes too eq_trajectoryFileName. Uses SuperCellSize to calculate volume."""
def eq_traj(myAtoms, trajObject, eq_trajObject, superCellSize):
    try:
        t = 0
        eq_index = 0
        while t < len(trajObject)-3: #Will check equilibrium conditions for atoms object and stop when equilibrium is reached.
            P_tot_diff = 0
            E_tot_diff = 0
            V_diff = 0
            for i in range(3): #Sums up offset of energy, pressure and volume between timesteps for three different values of i.
                curr_element = t+i
                next_element = t+i+1
                E_tot_1 = trajObject[curr_element].get_total_energy() #Total energy of two consecutive timesteps
                E_tot_2 = trajObject[next_element].get_total_energy()
                P_inst_1 = calc_instantaneous_pressure(myAtoms, trajObject, superCellSize, curr_element) #instantaneous pressure at two consecutive timesteps
                P_inst_2 = calc_instantaneous_pressure(myAtoms, trajObject, superCellSize, next_element) 
                V_1 = trajObject[curr_element].get_volume() * superCellSize #Volume of the cell at time two consecutive timesteps
                V_2 = trajObject[next_element].get_volume() * superCellSize 
                V_diff += abs(V_2 - V_1) #Offset in volume between timesteps.
                P_tot_diff += abs(P_inst_2 - P_inst_1) #Offset in pressure
                E_tot_diff += abs(E_tot_2 - E_tot_1) #Offset in total energy
            P_tot_diff_mean = P_tot_diff/3 #Mean values of three iterations
            V_diff_mean = V_diff/3
            E_tot_diff_mean = E_tot_diff/3
            if E_tot_diff_mean < 0.02 and V_diff_mean < 0.1 and P_tot_diff_mean < 1e-6 : #Criteria for equilibrium. Still not checking P_tot_diff_mean
                eq_index = t #saves index of first atom that has reached equilibrium.
                break
            t += 1
        t = len(trajObject) - 1
        for i in range(eq_index, len(trajObject)): #while loop that goes through all atoms objects in equilibrium and writes them to new .traj-file
            eq_trajObject.write(trajObject[i])
    except Exception as e:
        print("An error occured when checking equilibrium conditions", e)
        return(None)

# Calculates the specific heat and returns a numpy.float64 with dimensions J/(K*Kg)
def Specific_Heat(myAtoms, trajObject):   
    try:
        myAtoms.get_masses() # Tries if the attribute exists, skips except if it does
    except AttributeError:
        print("You have not entered a valid system.") # Message for user if no attribute
        return False # Ends the function
      
    bulk_mass=sum(myAtoms.get_masses())*1.6605402*10**(-27)
    temp_diff = (trajObject[1].get_kinetic_energy() /len(trajObject[1]) - trajObject[0].get_kinetic_energy() /len(trajObject[0])) / (1.5 * units.kB)  #Temperature difference between two runs when system has reached equilibrium
    pot_energy_diff = (trajObject[1].get_potential_energy() /len(trajObject[1]) 
                        - trajObject[0].get_potential_energy() /len(trajObject[0])) # potential energy difference when ystem has reached equilibrium

    kin_energy_diff = (trajObject[1].get_kinetic_energy() /len(trajObject[1]) 
                            - trajObject[0].get_kinetic_energy()/len(trajObject[0])) # potential energy difference when ystem has reached equilibrium
    
    heat_capcity = abs(((pot_energy_diff + kin_energy_diff)*(1.6021765*10**(-19)))/(temp_diff) / bulk_mass)
    print("C_p = ", heat_capcity, "[J/K*Kg]")
    return heat_capcity


"""Function to calculate and print the time average of the mean square displacement (MSD) at time t."""
def MSD_calc(myAtoms, trajObject, timeStepIndex):
    try:
        time = len(trajObject)-timeStepIndex
        pos_eq = trajObject[-1].get_positions() #position of atoms when system has reached equilibrium
        pos_t = trajObject[1].get_positions() #position of atoms at time t
        diff = pos_t - pos_eq #displacement of all atoms from equilibrium to time t as a vector
        diff_sq = np.absolute(diff)**2 
        MSD = np.sum(diff_sq)/len(myAtoms) #Time averaged mean square displacement.
    except Exception as e:
        print("An error occured when calculating the mean square displacement.", e)
        return(None)
    print("MSD = ", MSD, "[Å²]")
    return(MSD)

"""Function that  calculates the self-diffusion coefficient (D) at time t, based on the value of the mean square displacement."""

def Self_diffuse(trajObject, MSD, timeStepIndex):
    try:
        D = 5*MSD/(6*len(trajObject)) #How to connect mean squre displacement to self-diffusion coefficient. Multiply by 5 because timestep is 5 fs.
    except Exception as e:
        print("An error ocurred when calculating the self diffusion coefficient.", e)
        return(None)
    print("D = ", D, "[Å²/fs]")
    return(D)
    
"""Function that checks the Lindemann criterion which determines if the system is melting or not."""
def Lindemann(trajObject, MSD, timeStepIndex):
    try:
        nblist = FullNeighborList(3.5, trajObject[-1]).get_neighbors(1, -1) #Returns 3 lists containing information about nearest neighbors. 3rd list is the square of the distance to the neighbors.
        d = np.sqrt(np.amin(nblist[2])) #distance to the nearest neighbor. Takes the minimum value of nblist.
        L = MSD/d #Lindemann criterion. Expect melting when L>0.1
    except Exception as e:
        print("An error occured when checking the Lindemann criterion.", e)
        return(None)
    if L > 0.1:
            print("MELTING!")
            return True 
    else:
        print("NOT MELTING!")
        return False

def calc_instantaneous_pressure(myAtoms, trajObject, superCellSize, timeStepIndex):
    """ Calculates instantaneous pressure at time timeStepIndex * deltaT, i.e. P(n*deltaT).
    superCellSize is the repetition of unit cell in each direction in 3D space."""
    try: 
        atomsAtEquilibrium = trajObject[timeStepIndex]                  # Retrieve atoms some amount of timesteps in when system has reached equilibrium. 
        eqPos = atomsAtEquilibrium.get_positions()                      # Atom 3D coords. Unit of length [Å] (I think).
        eqTemperature = atomsAtEquilibrium.get_temperature()            # Temperature at time t = timeStepIndex * deltaT. Unit [K] (I think).
        eqVolume = atomsAtEquilibrium.get_volume() * superCellSize      # Get volume from atoms object at time t = timeStepIndex * deltaT. 
                                                                        # By getting the unit cell volume and scale with repetition. Unit [Å^3] I think.
                                                                        # This might be constant always but not sure for defect systems... 
        
        N = len(atomsAtEquilibrium)                                     # Number of atoms in object
        eqForces = myAtoms.calc.get_forces(atomsAtEquilibrium)          # Use calculator object to get the forces on the atoms. Unit: [eV/Å] (I think)
        
        # Instantaneos pressure 
        posForce = np.sum(eqPos * eqForces, axis=1)                     # Dot product of all forces and positions. Unit [eV * Å / Å]  = [eV]
        instantP = (N * units.kB * eqTemperature) / eqVolume + np.sum(posForce) / (3 * eqVolume) # Unit kB: [eV/K], -> [(eV/K) * K / Å^3 + eV / Å^3] = [eV / Å^3]
    except Exception as e:
        print("An error occured when calculating instantaneos pressure:", e)
        return(None)

    #print(instantP)
    return(instantP)

def calc_internal_pressure(myAtoms, trajObject, superCellSize):
    """ Internal pressure is the MD average of the instantaneous pressures. 
    IMPORTANT! trajObject must contain the atoms objects when the system has reached equilibrium,
    in order to get internal temperature from a stable crystal. 
    superCellSize is the repetition of unit cell in each direction in 3D space."""
    try:  
        iterationsInEqulibrium = len(trajObject)                        # Will be the amount of iterations in equilibrium
        M = iterationsInEqulibrium                                      # M = (iterations * deltaT) / deltaT
        
        # MD average
        allInstantPressures = 0
        for n in range(1,iterationsInEqulibrium):                       # Start from 1 since it is the first element in a .traj object. 
            instPressure = calc_instantaneous_pressure(myAtoms, trajObject, superCellSize,  n)
            if instPressure is not None:
                allInstantPressures += instPressure
            else:
                raise Exception("Instantenous pressure returned: None")
        
        internalPressure = allInstantPressures / M                      # Internal pressure is the MD average of the instantaneous pressures. 
                                                                        # Unit summation of instantaneosPressures => [eV / Å^3]
        print("Internal Pressure:", internalPressure, "[eV / Å^3]")
        return(internalPressure)
    except Exception as e:
        print("An error occured in internal pressure function:", e)
        return(None)
      
def calc_internal_temperature(myAtoms, trajectoryFileName, timeStepIndex):
    """ Returns the average temperature within parameters """
    
    for i in range(1,timeStepIndex):
        sumTemp = sum(myAtoms.get_temperature())     
        eqTemperature += sumTemp/len(sumTemp)               # Iteration gives the sum of a number of average temperatures

    avgTemperature = eqTemperature/timeStepIndex            # Sampling this gives the average temperature for the atoms
    print("Average internal temperature:", avgTemperature)  
    return(avgTemperature)

def calc_lattice_constant_fcc_cubic(atoms, atomName):
    """ Calculates the lattice constants. IMPORTANT!: Only works for FCC cubic crystals. Calculates both a and c constant but those are equal for fcc cubic"""
    # Make a good initial guess on the lattice constant
    a0 = 3.52 / np.sqrt(2) 
    c0 = np.sqrt(8 / 3.0) * a0
    
    fileName = "lattice_" + atomName + ".traj"                              # Create filename from atomname
    traj = Trajectory(fileName, 'w')                                        # Create a traj file to store the results from calculations

    # Generate 9 calculations of potential energy for different a and c values. 
    eps = 0.01                                                              # A small deviation to generate a few more constants
    for a in a0 * np.linspace(1 - eps, 1 + eps, 3):
        for c in c0 * np.linspace(1 - eps, 1 + eps, 3):
            at = bulk(atomName, 'fcc', a=a, c=c, cubic=True)                # Use bulk to build a cell. Only config is cubic fcc.
            #at.calc = atoms.calc                                           # Something doesn't work with copying calculator from atoms object. 
            at.calc = EMT()                                                 # Assign a new calculator.
            at.get_potential_energy()                                       
            traj.write(at)                                                  # Write both bulk config and energy to trajectory file

    # Now we can get the energies and lattice constants from the traj file
    configs = read(fileName + "@:")
    energies = [config.get_potential_energy() for config in configs]        # Get the atoms objects from traj file. 

    # From the bulk builder we can do at.cell to get the constant: a at position 0,0 and constant: c at position 2,2
    a = np.array([config.cell[0, 0] for config in configs])                 
    c = np.array([config.cell[2, 2] for config in configs])

    # Fit the energy to lattice constants to the expression: 𝑝0+𝑝1𝑎+𝑝2𝑐+𝑝3𝑎^2+𝑝4𝑎𝑐+𝑝5𝑐^2
    functions = np.array([a**0, a, c, a**2, a * c, c**2])
    p = np.linalg.lstsq(functions.T, energies, rcond=-1)[0]
    #print("Polynomial:", p, "\n")
    #print("Functions: \n", functions.T)
    
    # Solve fitted function for a and c. The minimum is found by
    p0 = p[0]
    p1 = p[1:3]
    p2 = np.array([(2 * p[3], p[4]), (p[4], 2 * p[5])])
    a0, c0 = np.linalg.solve(p2.T, -p1)

    print("Lattice constant a:", a0) 
    return(a0)
    #print("Lattice constants a:", a0, "| c:", c0, "\n") Uncomment if we want to print c also
