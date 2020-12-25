""" Initiation of variables and system  """
import math
import shutil
# 6 of the 7 lattice systems (rhombohedral is not available)
import ase.io
import ase
from ase.lattice.cubic import *
from ase.lattice.tetragonal import *
from ase.lattice.orthorhombic import *
from ase.lattice.monoclinic import *
from ase.lattice.triclinic import *
from ase.lattice.hexagonal import *
from ase.atom import *
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

# Dependencies to run materials project
from pymatgen.ext.matproj import MPRester
from pymatgen.io.cif import CifParser

# Initiation functions to separate them from variables
from .init_functions import create_vacancy, find_crystal_center, set_lattice
from .init_functions import set_lattice_const, insert_impurity

atoms_list = []
atoms_opt_list = []

# Check if the potential string is empty, then assign standard Lennard Jones potential
def checkKIMpotential(potential):
    if (potential == "") or (potential == " "):
        return "LJ_ElliottAkerson_2015_Universal__MO_959249795837_003"
    else:
        return potential

# Optimized lattice const initiation

def init_opt(EMT_Check, openKIM_Check, Lennard_Jones_Check, LJ_epsilon,
                            LJ_sigma, LJ_cutoff,Verlocity_Verlet_Check, KIM_potential,Symbol,
                            Vacancy, Impurity, Impurity_ele, Temperature,
                            Size_X,Size_Y,Size_Z,PBC,Bravais_lattice,Directions,Miller,
                            lc_a,lc_b,lc_c,lc_alpha,lc_beta,lc_gamma):

    Lattice_Const = set_lattice_const(lc_a,
                                    lc_b,
                                    lc_c,
                                    lc_alpha,
                                    lc_beta,
                                    lc_gamma)

    # Set up a crystals
    atoms_opt = set_lattice(getattr(ase.lattice.cubic, Bravais_lattice),
                    Lattice_Const,
                    Directions,
                    Miller,
                    Size_X,
                    Size_Y,
                    Size_Z,
                    Symbol,
                    PBC) 

    #places impurity in the crystal 
    if Impurity == True:
        atom_pos = find_crystal_center(atoms_opt) # Returns a center position in the crystal
        insert_impurity(atoms_opt, Impurity_ele, atom_pos) # Insert "foregin" atom in the crystal

    #Places vacancy in the crytal
    if Vacancy == True:
        create_vacancy(atoms_opt) # Create a vacancy

    # Set the momenta corresponding to desired temperature when running Verlocity Verlet
    if Verlocity_Verlet_Check == True:
        MaxwellBoltzmannDistribution(atoms_opt, Temperature * units.kB)
        Stationary(atoms_opt) # Set linear momentum to zero
        ZeroRotation(atoms_opt) # Set angular momentum to zero

    # Interatomic potential
    if (EMT_Check == True) and (openKIM_Check == False) and (Lennard_Jones_Check == False):
        atoms_opt.calc = EMT()
    elif (EMT_Check == False) and (openKIM_Check == True) and (Lennard_Jones_Check == False):
        #Sets the potential for openKIM. If none is given returns standard Lennard-Jones
        potential = checkKIMpotential(KIM_potential)
        atoms_opt.calc = KIM(potential, options={"ase_neigh": True})
        #atoms.set_calculator(OpenKIMcalculator(potential))
    elif (EMT_Check == False) and (openKIM_Check == False) and (Lennard_Jones_Check == True):
        atoms_opt.calc = LennardJones(list(dict.fromkeys(atoms.get_atomic_numbers())), LJ_epsilon, LJ_sigma, rCut=LJ_cutoff, modified=True)
    else:
        raise Exception("Only one of EMT, OpenKim and Lennard_jones can be true!")

    atoms_opt_list.append(atoms_opt)
    return atoms_opt_list

# Init for ASE
def init(EMT_Check, openKIM_Check, Lennard_Jones_Check, LJ_epsilon,
                            LJ_sigma, LJ_cutoff,Verlocity_Verlet_Check, KIM_potential,Symbol,
                            Vacancy, Impurity, Impurity_ele, Temperature,
                            Size_X,Size_Y,Size_Z,PBC,Bravais_lattice,Directions,Miller,
                            lc_a,lc_b,lc_c,lc_alpha,lc_beta,lc_gamma):

    Lattice_Const = set_lattice_const(lc_a,
                                    lc_b,
                                    lc_c,
                                    lc_alpha,
                                    lc_beta,
                                    lc_gamma)

    # Set up a crystals
    atoms = set_lattice(getattr(ase.lattice.cubic, Bravais_lattice),
                    Lattice_Const,
                    Directions,
                    Miller,
                    Size_X,
                    Size_Y,
                    Size_Z,
                    Symbol,
                    PBC) 

    #places impurity in the crystal 
    if Impurity == True:
        atom_pos = find_crystal_center(atoms) # Returns a center position in the crystal
        insert_impurity(atoms, Impurity_ele, atom_pos) # Insert "foregin" atom in the crystal

    #Places vacancy in the crytal
    if Vacancy == True:
        create_vacancy(atoms) # Create a vacancy

    # Set the momenta corresponding to desired temperature when running Verlocity Verlet
    if Verlocity_Verlet_Check == True:
        MaxwellBoltzmannDistribution(atoms, Temperature * units.kB)
        Stationary(atoms) # Set linear momentum to zero
        ZeroRotation(atoms) # Set angular momentum to zero

    # Interatomic potential
    if (EMT_Check == True) and (openKIM_Check == False) and (Lennard_Jones_Check == False):
        atoms.calc = EMT()
    elif (EMT_Check == False) and (openKIM_Check == True) and (Lennard_Jones_Check == False):
        #Sets the potential for openKIM. If none is given returns standard Lennard-Jones
        potential = checkKIMpotential(KIM_potential)
        atoms.calc = KIM(potential, options={"ase_neigh": True})
        #atoms.set_calculator(OpenKIMcalculator(potential))
    elif (EMT_Check == False) and (openKIM_Check == False) and (Lennard_Jones_Check == True):
        atoms.calc = LennardJones(list(dict.fromkeys(atoms.get_atomic_numbers())), LJ_epsilon, LJ_sigma, rCut=LJ_cutoff, modified=True)
    else:
        raise Exception("Only one of EMT, OpenKim and Lennard_jones can be true!")

    atoms_list.append(atoms)
    return atoms_list

# Init for Materials project
def init_MP(EMT_Check,openKIM_Check,Lennard_Jones_Check, LJ_epsilon,
                                LJ_sigma, LJ_cutoff,Verlocity_Verlet_Check,KIM_potential,Criteria_list,
                                Vacancy, Impurity, Impurity_ele, Temperature,
                                Size_X,Size_Y,Size_Z,API_Key,PBC):
    m = MPRester(API_Key) 
    
    #Loop that takes out each critera for each query
    for criteria in Criteria_list:

        #If there are no elements in data raise an exception and end program
        data = m.query(criteria, properties=['cif', 'spacegroup', 'pretty_formula'])
        print(data)
        if len(data) != 0:
            for i in range(len(data)):
            
                # Takes out the space group and crystal structure query
                space_group = ((data[i])['spacegroup'])['symbol']
                crystal_structure = ((data[i])['spacegroup'])['crystal_system']
                pretty_formula = str((data[i])['pretty_formula'])
            
                # Function that skips the element if it not an FCC crystal
                if space_group[0] != 'F' or crystal_structure != 'cubic':
                    continue
                
                #Takes out the CIF information and creates a unique file
                cif_Info=(data[i])["cif"]
                f = open(pretty_formula + ".cif", "w+")
                f.write(cif_Info)
                f.close()
                
                #Creates the atom object from the CIF information
                atoms = ase.io.read(pretty_formula + ".cif")
                atoms.set_pbc(PBC)
                
                #Creates a supercell
                atoms = atoms*(Size_X,Size_Y,Size_Z)

                #places impurity in the crystal 
                if Impurity == True:
                    atom_pos = find_crystal_center(atoms) # Returns a center position in the crystal
                    insert_impurity(atoms, Impurity_ele, atom_pos) # Insert "foregin" atom in the crystal

                #Places vacancy in the crytal
                if Vacancy == True:
                    create_vacancy(atoms) # Create a vacancy

                # Set the momenta corresponding to desired temperature when running Verlocity Verlet
                if Verlocity_Verlet_Check == True:
                    MaxwellBoltzmannDistribution(atoms, Temperature * units.kB)
                    Stationary(atoms) # Set linear momentum to zero
                    ZeroRotation(atoms) # Set angular momentum to zero
                # Interatomic potential
                if (EMT_Check == True) and (openKIM_Check == False) and (Lennard_Jones_Check == False):
                    atoms.calc = EMT()
                elif (EMT_Check == False) and (openKIM_Check == True) and (Lennard_Jones_Check == False):
                    #Sets the potential for openKIM. If none is given returns standard Lennard-Jones
                    potential = checkKIMpotential(KIM_potential)
                    atoms.calc = KIM(potential, options={"ase_neigh": True})
                    #atoms.set_calculator(OpenKIMcalculator(potential))
                elif (EMT_Check == False) and (openKIM_Check == False) and (Lennard_Jones_Check == True):
                        atoms.calc = LennardJones(list(dict.fromkeys(atoms.get_atomic_numbers())), LJ_epsilon, LJ_sigma, rCut=LJ_cutoff, modified=True)

                else:
                    raise Exception("Only one of EMT, OpenKim and Lennard_jones can be true!")
                
                #Moves the trajectory file to another folder after it has been used
                shutil.move(pretty_formula + ".cif", "CIF/" + pretty_formula + ".cif")
                atoms_list.append(atoms)
        else:
            raise Exception("The query/queries returned no elements!") 
    
    # If atoms_list is empty raise an expetion and end program
    if len(atoms_list) == 0:
        raise Exception("The query/queries contains no FCC crystals")
    else:
        return atoms_list
