from ase import units
from ase.data import atomic_masses, atomic_numbers
import numpy as np
from asap3 import Trajectory, FullNeighborList

from ase.calculators.eam import EAM
from ase.build import bulk
from ase.io import read
import sys, os
from ase.units import kJ
from ase.eos import EquationOfState

import csv
import datetime

"""Function that takes all the atoms-objects after the system reaches equilibrium  (constant total energy, volume and pressure) and writes them over to a new .traj-file. Goes through trajectoryFileName and writes too eq_trajectoryFileName. Uses SuperCellSize to calculate volume."""
def eq_traj(myAtoms, trajObject, superCellSize):
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
            print(E_tot_diff_mean)
            if E_tot_diff_mean < 0.001 : #Criteria for equilibrium. Still not checking P_tot_diff_mean
                eq_index = t #saves index of first atom that has reached equilibrium.
                break
            t += 1
        return eq_index
    except Exception as e:
        print("An error occured when calculating the checking the conditions for equilibrium:")
        exc_type, exc_obj, exc_traceBack = sys.exc_info()
        fname = os.path.split(exc_traceBack.tb_frame.f_code.co_filename)[1]
        print("Error type:", exc_type, "; Message:", e, "; In file:", fname, "; On line:", exc_traceBack.tb_lineno)
        return(None)

# Calculates the specific heat and returns a numpy.float64 with dimensions J/(K*Kg)
def Specific_Heat(myAtoms, trajObject, eq_index):   
    try:
        myAtoms.get_masses() # Tries if the attribute exists, skips except if it does
    except AttributeError:
        print("You have not entered a valid system.") # Message for user if no attribute
        return False # Ends the function
    eq_next = eq_index + 1  
    bulk_mass=sum(myAtoms.get_masses())*units._amu
    temp_diff = (trajObject[eq_next].get_kinetic_energy() /len(trajObject[eq_next]) - trajObject[eq_index].get_kinetic_energy() /len(trajObject[eq_index])) / (1.5 * units.kB)  #Temperature difference between two runs when system has reached equilibrium
    pot_energy_diff = (trajObject[eq_next].get_potential_energy() /len(trajObject[eq_next]) 
                        - trajObject[eq_index].get_potential_energy() /len(trajObject[eq_index])) # potential energy difference when ystem has reached equilibrium

    kin_energy_diff = (trajObject[eq_next].get_kinetic_energy() /len(trajObject[eq_next]) 
                            - trajObject[eq_index].get_kinetic_energy()/len(trajObject[eq_index])) # potential energy difference when ystem has reached equilibrium
    
    heat_capacity = abs(((pot_energy_diff + kin_energy_diff)*(1.6021765*10**(-19)))/(temp_diff) / bulk_mass)
    return heat_capacity


"""Function to calculate and print the time average of the mean square displacement (MSD) at time t."""
def MSD_calc(myAtoms, trajObject, timeStepIndex):
    try:
        pos_t = trajObject[-1].get_positions() #position of atoms when system has reached equilibrium
        pos_eq = trajObject[timeStepIndex].get_positions() #position of atoms at time t
        diff = pos_t - pos_eq #displacement of all atoms from equilibrium to time t as a vector
        diff_sq = np.absolute(diff)**2 
        MSD = np.sum(diff_sq)/len(myAtoms) #Time averaged mean square displacement.
    except Exception as e:
        print("An error occured when calculating the mean square displacement:")
        exc_type, exc_obj, exc_traceBack = sys.exc_info()
        fname = os.path.split(exc_traceBack.tb_frame.f_code.co_filename)[1]
        print("Error type:", exc_type, "; Message:", e, "; In file:", fname, "; On line:", exc_traceBack.tb_lineno)
        return(None)
    return(MSD)

"""Function that  calculates the self-diffusion coefficient (D) at time t, based on the value of the mean square displacement."""

def Self_diffuse(MSD, timeStepIndex):
    try:
        D = 5*MSD/(6*timeStepIndex) #How to connect mean squre displacement to self-diffusion coefficient. Multiply by 5 because timestep is 5 fs.
    except Exception as e:
        print("An error occured when calculating the self diffusion coefficient:")
        exc_type, exc_obj, exc_traceBack = sys.exc_info()
        fname = os.path.split(exc_traceBack.tb_frame.f_code.co_filename)[1]
        print("Error type:", exc_type, "; Message:", e, "; In file:", fname, "; On line:", exc_traceBack.tb_lineno)
        return(None)
    return(D)
    
"""Function that checks the Lindemann criterion which determines if the system is melting or not."""
def Lindemann(trajObject, MSD):
    try:
        nblist = FullNeighborList(3.5, trajObject[-1]).get_neighbors(1, -1) #Returns 3 lists containing information about nearest neighbors. 3rd list is the square of the distance to the neighbors.
        d = np.sqrt(np.amin(nblist[2])) #distance to the nearest neighbor. Takes the minimum value of nblist.
        L = MSD/d #Lindemann criterion. Expect melting when L>0.1
    except Exception as e:
        print("An error occured when checking the Lindemann criterion:")
        exc_type, exc_obj, exc_traceBack = sys.exc_info()
        fname = os.path.split(exc_traceBack.tb_frame.f_code.co_filename)[1]
        print("Error type:", exc_type, "; Message:", e, "; In file:", fname, "; On line:", exc_traceBack.tb_lineno)
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

def calc_internal_pressure(myAtoms, trajObject, eq_index, superCellSize):
    """ Internal pressure is the MD average of the instantaneous pressures. 
    IMPORTANT! trajObject must contain the atoms objects when the system has reached equilibrium,
    in order to get internal temperature from a stable crystal. 
    superCellSize is the repetition of unit cell in each direction in 3D space."""
    try:  
        iterationsInEqulibrium = len(trajObject) - eq_index                        # Will be the amount of iterations in equilibrium
        M = iterationsInEqulibrium                                      # M = (iterations * deltaT) / deltaT
        
        # MD average
        allInstantPressures = 0
        for n in range(eq_index, len(trajObject)):                       # Start from 1 since it is the first element in a .traj object. 
            instPressure = calc_instantaneous_pressure(myAtoms, trajObject, superCellSize, n)
            if instPressure is not None:
                allInstantPressures += instPressure
            else:
                raise Exception("Instantenous pressure returned: None")
        
        internalPressure = allInstantPressures / M                      # Internal pressure is the MD average of the instantaneous pressures. 
                                                                        # Unit summation of instantaneosPressures => [eV / Å^3]
        return(internalPressure)
    except Exception as e:
        print("An error occured in internal pressure function:", e)
        return(None)

def internal_temperature(atoms, trajObject, eq_index):
    """
    Calculates the internal temperature of the system
    
    Parameters
        atoms       :   Lattice of atoms being simulated
        trajObject  :   TrajectoryReader

    Returns the average of a sum of samples using the equation T = 2K/3kB,
    where K is the kinetic energy per atom
    """

    try:
        eq_length = len(trajObject) - eq_index
        eqEkin = 0
        for n in range(eq_index, len(trajObject)):                       
            eqEkin += trajObject[n].get_kinetic_energy()/len(atoms)                          # Sum kinetic energies for each trajectory object
        avgTemp = (2*eqEkin)/(3*units.kB*eq_length)                          # Average sum over number of samples and calculate temperature

    except Exception as e:
        print("An error occured when calculating the internal temperature")
        exc_type, exc_obj, exc_traceBack = sys.exc_info()
        fname = os.path.split(exc_traceBack.tb_frame.f_code.co_filename)[1]
        print("Error type:", exc_type, "; Message:", e, "; In file:", fname, "; On line:", exc_traceBack.tb_lineno)
        return(None)
    
    return(avgTemp)
    
def cohesive_energy(atoms, trajObject, eq_index):
    """
    Returns the cohesive energy of the system

    Parameters
        atoms       :   Lattice of atoms being simulated
        trajObject  :   TrajectoryReader

    Returns the average of a sum of samples over the potential energy per atom
    """

    try:
        eq_length = len(trajObject) - eq_index
        eqEcoh = 0
        for n in range(eq_index, len(trajObject)):
            eqEcoh += trajObject[n].get_potential_energy()/len(atoms)           # Sum potential energies per atom for each trajectory object 
        avgEcoh = eqEcoh/eq_length                                        # Average sum over number of samples
    
    except Exception as e:
        print("An error occured when calculating the cohesive energy")
        exc_type, exc_obj, exc_traceBack = sys.exc_info()
        fname = os.path.split(exc_traceBack.tb_frame.f_code.co_filename)[1]
        print("Error type:", exc_type, "; Message:", e, "; In file:", fname, "; On line:", exc_traceBack.tb_lineno)
        return(None)
    
    return(avgEcoh)

def debye_temperature(trajObject, MSD):
    """
    Calculates the Debye temperature of the system.

    Parameters
        trajObject  : TrajectoryReader
        MSD         : Returned from MSD_calc()

    Uses the functions for temperature, atom masses converted to kg and mean square displacement 
    Returns the average of a sum of samples over the Debye temperature of the system
    """
    try: 
        eqDebye = 0
        for n in range(1, len(trajObject)):
            T = trajObject[n].get_temperature()                                # Set to system temperature                     
            m = sum(trajObject[n].get_masses())*units._amu                     # Set to sum of atom masses converted to kg
            eqDebye += np.sqrt((3*(units._hbar**2)*T)/(m*units.kB*MSD))        # Sum Debye temperatures for each trajectory object
        avgDebye = eqDebye/len(trajObject)                                     # Average sum over number of samples
    
    except Exception as e:
        print("An error occured when calculating the Debye temperature:")
        exc_type, exc_obj, exc_traceBack = sys.exc_info()
        fname = os.path.split(exc_traceBack.tb_frame.f_code.co_filename)[1]
        print("Error type:", exc_type, "; Message:", e, "; In file:", fname, "; On line:", exc_traceBack.tb_lineno)
        return(None)

    return(avgDebye)

def calc_lattice_constant_fcc_cubic(atomName, atomsCalculator):
    """ Calculates the lattice constants. IMPORTANT!: Only works for FCC cubic crystals. 
        Calculates both a and c constant but those are equal for fcc cubic.
        Only calculates for pure one atom crystals. Modification for defect systems might have to be made, using the original atoms object maybe."""
    try: 
        # Make a good initial guess on the lattice constant
        a0 = 3.52 / np.sqrt(2) 
        c0 = np.sqrt(8 / 3.0) * a0
        
        fileName = "lattice_" + atomName + ".traj"                              # Create filename from atomname.
        traj = Trajectory(fileName, 'w')                                        # Create a traj file to store the results from calculations.

        # Generate 9 calculations of potential energy for different a and c values. 
        eps = 0.01                                                              # A small deviation to generate a few more constants.
        for a in a0 * np.linspace(1 - eps, 1 + eps, 3):
            for c in c0 * np.linspace(1 - eps, 1 + eps, 3):
                at = bulk(atomName, 'fcc', a=a, c=c, cubic=True)                # Use bulk to build a cell. Only config is cubic fcc.
                at.calc = atomsCalculator                                       # Assign calculator that is in original atoms object.                        
                traj.write(at)                                                  # Write bulk config to trajectory file

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

        return(a0)
        #print("Lattice constants a:", a0, "| c:", c0, "\n") Uncomment if we want to print c also
    except Exception as e:
        print("An error occured when calculating the lattice constant:")
        exc_type, exc_obj, exc_traceBack = sys.exc_info()
        fname = os.path.split(exc_traceBack.tb_frame.f_code.co_filename)[1]
        print("Error type:", exc_type, "; Message:", e, "; In file:", fname, "; On line:", exc_traceBack.tb_lineno)
        return(None)

def calc_bulk_modulus(atoms):
    """
    Calculates the bulk modulus using Equation of State (EOS). Based on the volume, this can also be used for calculation of lattice constant for cubic crystal structures. 
    """
    try:
        atomsConfigName = str(atoms.symbols)                    # Convert symbolsname to string, will make it as a 'molecule' notation, i.e chemical symbols + amount of atoms. 
        trajFileName = "Bulk_" + atomsConfigName + '.traj'
        cell = atoms.get_cell()                                 # Extract atoms cell, i.e. lattice constants to make small deviations in for-loop.
        traj = Trajectory(trajFileName, 'w')                    # Create a traj file for writing results to.
        measAmount = 10                                         # This is set to 10 from experimenting with different materials that yields an energy minimum between endpoints. 
                                                                # Please note that this can be modified to a larger value if the curve from EOS is a straight line for instance.
        for x in np.linspace(0.8, 1.2, measAmount):             # Generates a array of measAmount numbers equally spaced between 0.8 and 1.2 to vary the lattice constants with
            atoms.set_cell(cell * x, scale_atoms=True)          # Modify the cell with new lattice parameters, to later plot energy to volume. Careful this modifies the original object, must reset below. 
            traj.write(atoms)
        
        atoms.set_cell(cell, scale_atoms=True)                  # Reset the cell to original. 

        configs = read(trajFileName + '@:')                     # To read all configs in traj file @ followed by a intervall. Example 0:4 first 5 elements, ':' denotes all elements. 
                                                                # measAmount elements from the for loop above
        
        # Extract volumes and energies, divide with amount of atoms to scale for a cell
        volumes = [atoms.get_volume()/len(atoms) for atoms in configs]                 
        energies = [atoms.get_potential_energy()/len(atoms) for atoms in configs]
        eos = EquationOfState(volumes, energies)                # Generate EOS
        v0, e0, B = eos.fit()                                   # This returns minimized energy e0 for the corresponding volume v0, and B is the bulk modulus (curvature at v0).
        B_GPa = B / kJ * 1.0e24                                 # Unit conversion

        if not (e0 < energies[0] and e0 < energies[-1]):        # Have to check that minmum is not an endpoint. Can replicate with bad lattice constant guess. 
            raise ValueError("Minumum is endpoint, use a different intervall. Or make a better guess on lattice constant.")
 
        eos.plot(atomsConfigName + '_eos.png')                  # Saves an images of the plot. Might not want this?

        return(e0, v0, B_GPa)                                   # Return min E=e0 at volume V=v0 and Bulkmodulus with unit GPa

    except Exception as e:
        print("An error occured when calculating the Bulk modulus:")
        exc_type, exc_obj, exc_traceBack = sys.exc_info()
        fname = os.path.split(exc_traceBack.tb_frame.f_code.co_filename)[1]
        print("Error type:", exc_type, "; Message:", e, "; In file:", fname, "; On line:", exc_traceBack.tb_lineno)
        return(None, None, None)

def write_atom_properties(myAtoms, csvFileName, trajObject, eq_index):
    """calculates a set of chosen properties and saves them to a csv-file (comma seperated values). The file is to be used to make plots of the results."""
    try:
        eq_length = len(trajObject) - eq_index
        file = open(csvFileName, "w", newline="")
        fieldnames = ["Time", "MSD", "S"] #These will be the first rows in each column specifying the property."
        writer = csv.DictWriter(file, fieldnames=fieldnames, delimiter=";") #makes a dictionary-writer for csv-files.
        t = 1
        writer.writeheader() #Sets fieldnames as first rows in file.
        while t < eq_length: #Checks time evolution of chosen properties.
            time_after_eq = eq_length - t
            MSD = MSD_calc(myAtoms, trajObject, time_after_eq)
            S = Self_diffuse(MSD, t)
            writer.writerow({"Time" : t, "MSD" : MSD, "S" : S}) #Writes values at time t to csv-file. A new row is a new timestep.
            t += 1
    except Exception as e:
        print("An error occured when writing values to .csv-file:")
        exc_type, exc_obj, exc_traceBack = sys.exc_info()
        fname = os.path.split(exc_traceBack.tb_frame.f_code.co_filename)[1]
        print("Error type:", exc_type, "; Message:", e, "; In file:", fname, "; On line:", exc_traceBack.tb_lineno)
