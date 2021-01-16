# Initiate the crystal based on the chosen variables
# This will eventually become "Initiate the system" => system depends on user's choice
import os
import shutil
from Init.init_values import *
from tkinter import *
import Json.write_json as write_json
from mpi4py import MPI

from Init.init_functions import create_vacancy, find_crystal_center, insert_impurity
from Init.init_values import checkKIMpotential
import Calculations.calculations_mpi as calc
from Optimade.optimade import translate_to_optimade, concatenateOptimadeDataFiles
from translate import translate

from ase.gui import *
from ase import units
from ase.calculators.kim.kim import KIM

# Algorithms and calculators for the simulation
from asap3.md.velocitydistribution import MaxwellBoltzmannDistribution, Stationary
from ase.md.velocitydistribution import ZeroRotation
from asap3.md.verlet import VelocityVerlet
from asap3.md.langevin import Langevin
from asap3 import EMT
from asap3 import Trajectory
from asap3 import LennardJones
#from asap3 import OpenKIMcalculator 

os.environ['OPENBLAS_NUM_THREADS'] = "1"
os.environ['MKL_NUM_THREADS'] = "1"
os.environ['NUMEXPR_NUM_THREADS'] = "1"
os.environ['OMP_NUM_THREADS'] = "1"

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

""" Ignore the parallel rank deprecation warning """

def warn(*args, **kwargs):
    pass
import warnings
warnings.warn = warn

def main():
    if rank == 0:
        print("We have", size, "processors.")

        EMT_Check,openKIM_Check, Lennard_Jones_Check, LJ_epsilon, LJ_sigma, LJ_cutoff, Velocity_Verlet_Check, Langevin_Check, Langevin_friction, Time_step, KIM_potential, ASE, Symbol, Materials_project,API_Key,Criteria_list, Vacancy, Impurity, Impurity_ele, Temperature, Steps, Interval,Size_X, Size_Y, Size_Z, PBC, Bravais_lattice,Directions,Miller, lc_a,lc_b,lc_c,lc_alpha,lc_beta,lc_gamma,Run_Optimade,Optimade_name,Optimized_volume = translate()

        """ Function that looks if the user wants to run ASE or Materials_project 
            Checks if the simulation is going to add impurites or not
            Impurites doesnt work with openKIM however
        """
        if (ASE == True) and (Materials_project == False):

            atoms = init(Symbol, Bravais_lattice, Directions, Miller,
                            lc_a,lc_b,lc_c,lc_alpha,lc_beta,lc_gamma)


        elif (Materials_project == True) and (ASE == False):
            atoms = init_MP(Criteria_list, API_Key)

        else:
            raise Exception("ASE=Materials_Project. Both cannot be true/false at the same time!")
        
        # atoms_array = np.array(atoms)
        # job_array = np.array_split(atoms_array, size)
        for i in range(0, size):
         #   print(job_array[i])
            comm.isend(atoms, dest=i, tag=i)
    
    my_jobs = comm.recv(source=0, tag=rank)
    print(my_jobs)

 #   result = [run_algorithm(atomobj[0]) for atomobj in my_jobs]
    result = run_algorithm(my_jobs[rank])
    comm.isend(result, dest=0, tag=0)
    if rank == 0:
        for i in range(0, size):
            result = comm.recv(source=i, tag=0)


def run_algorithm(atomobj):

        # # Run simulation with optimized volume.
        # if Optimized_volume == True:

        #     # Save all the angles for the unit cell
        #     alpha = (atomobj.get_cell_lengths_and_angles())[3]
        #     beta = (atomobj.get_cell_lengths_and_angles())[4]
        #     gamma = (atomobj.get_cell_lengths_and_angles())[5]

        #     #Calculates the optimal length for the lattice constant 
        #     if EMT_Check == True:
        #         latticeConstant_c = calc.calc_lattice_constant_cubic(atomobj, EMT(), alpha, beta, gamma, Size_X, Size_Y, Size_Z, PBC)
        #     elif openKIM_Check == True:
        #         potential = checkKIMpotential(KIM_potential)
        #         latticeConstant_c = calc.calc_lattice_constant_cubic(atomobj, KIM(potential, options={"ase_neigh": True}), 
        #                                                                                  alpha, beta, gamma, Size_X, Size_Y, Size_Z, PBC)
        #         # latticeConstant_c = calc.calc_lattice_constant_cubic(atomobj, OpenKIMcalculator(potential), 
        #         #                                                                          alpha, beta, gamma, Size_X, Size_Y, Size_Z, PBC)
        #     elif Lennard_Jones_Check == True:
        #         latticeConstant_c = calc.calc_lattice_constant_cubic(atomobj, LennardJones(list(dict.fromkeys(atomobj.get_atomic_numbers())), 
        #                                                                                 LJ_epsilon, LJ_sigma, rCut=LJ_cutoff, modified=True), 
        #                                                                                 alpha, beta, gamma, Size_X, Size_Y, Size_Z, PBC)
        #     # Updates the cell with the optimal lattice constant
        #     atomobj.set_cell([latticeConstant_c, latticeConstant_c, latticeConstant_c, alpha, beta, gamma])
        #     print("lattice constant:", latticeConstant_c, "\n")
        # else:
        #     # The lattice constant is set to zero if Optimized_volume == false. This is implemented so that the visualization does not crash!
        #     latticeConstant_c = 0 

        #Creates a supercell
    atomobj = atomobj*(5,5,5)

        #Set the periodic boundary conditions
    atomobj.set_pbc([True,True,True])

        # #places impurity in the crystal 
        # if Impurity == True:
        #     atom_pos = find_crystal_center(atomobj) # Returns a center position in the crystal
        #     insert_impurity(atomobj, Impurity_ele, atom_pos) # Insert "foregin" atom in the crystal

        # #Places vacancy in the crytal
        # if Vacancy == True:
        #     create_vacancy(atomobj) # Create a vacancy

        # Set the momenta corresponding to desired temperature when running Verlocity Verlet
        # if Velocity_Verlet_Check == True:
    MaxwellBoltzmannDistribution(atomobj, 300 * units.kB)
    Stationary(atomobj) # Set linear momentum to zero
    ZeroRotation(atomobj) # Set angular momentum to zero
        # # Interatomic potential
        # if (EMT_Check == True) and (openKIM_Check == False) and (Lennard_Jones_Check == False):
    atomobj.calc = EMT()
        # elif (EMT_Check == False) and (openKIM_Check == True) and (Lennard_Jones_Check == False):
        #     #Sets the potential for openKIM. If none is given returns standard Lennard-Jones
        #     potential = checkKIMpotential(KIM_potential)
        #     atomobj.calc = KIM(potential, options={"ase_neigh": True})
        #     #atomobj.set_calculator(OpenKIMcalculator(potential))
        # elif (EMT_Check == False) and (openKIM_Check == False) and (Lennard_Jones_Check == True):
        #     atomobj.calc = LennardJones(list(dict.fromkeys(atomobj.get_atomic_numbers())), LJ_epsilon, LJ_sigma, rCut=LJ_cutoff, modified=True)

        # else:
        #     raise Exception("Only one of EMT, OpenKim and Lennard_jones can be true!")

        # if (Velocity_Verlet_Check == True) and (Langevin_Check == False):

            # We want to run MD with constant energy using the VelocityVerlet algorithm.
    dyn = VelocityVerlet(atomobj, 1000*units.fs)
            
        # elif (Velocity_Verlet_Check == False) and (Langevin_Check == True):
        #     dyn = Langevin(atomobj, Time_step*units.fs, units.kB*Temperature, Langevin_friction)
        # else:
        #     raise Exception("Velocity_Verlet=Langevin. Both cannot be true/false at the same time!")
        #Creates a unique name for every simulation run 

        #Averaged kinetic energy, kinetic energy squared and temperature calculated
    k = []
    k_squared = []
    temp = []
    Steps = 1000
    Interval = 10
    
    for i in range(int(Steps/Interval)):
        dyn.run(Interval)
        k.append(atomobj.get_kinetic_energy())
        k_squared.append(atomobj.get_kinetic_energy()**2)
        temp.append(atomobj.get_temperature()) 
    avg_k_squared = sum(k_squared)/len(k_squared)
    avg_k = sum(k)/len(k)
    avg_temp = sum(temp)/len(temp)

    Nr_of_atoms = len(atomobj)

    bulk_mass=sum(atomobj.get_masses())*units._amu

    # if Velocity_Verlet_Check == True:    
    specificHeatCapacity = calc.Heat_Capcity_NVE(avg_k, avg_k_squared, avg_temp, bulk_mass, Nr_of_atoms)
    # else:
    #     specificHeatCapacity = calc.Heat_Capcity_NVT(atomobj)

    print("C_v = ", specificHeatCapacity, "[J/K*Kg]")
    print("rank" + str(rank) + "=done")



if __name__ == "__main__":
    main()