""" Initiation functions """

from collections import Counter
from ase import Atoms

# Parameters will decide values and bravais lattice
def set_lattice(Bravais,
                Lattice_Const,
                Directions,
                Miller,
                Size_X,
                Size_Y,
                Size_Z,
                Symbol,
                Pbc):
    
    atoms = Bravais(latticeconstant=Lattice_Const,
                    directions=Directions,
                    miller=Miller,
                    symbol=Symbol,
                    size=(Size_X, Size_Y, Size_Z),
                    pbc=Pbc)
    return atoms

# Handle choice of lattice constants
def set_lattice_const(lc_a, lc_b, lc_c, lc_alpha, lc_beta, lc_gamma):
    
    # Make a collection
    lc_constants = Counter(a = lc_a, 
                           b = lc_b, 
                           c = lc_c, 
                           alpha = lc_alpha,
                           beta = lc_beta,
                           gamma = lc_gamma)
    
    # Get the dictionary and remove zeroes
    lc_constants = dict(lc_constants)
    lc_constants = {x:y for x,y in lc_constants.items() if y!=0}
    
    # Checks if only lc_a is given since Bravais Cubic methods don't accept dict-format
    if (len(lc_constants) == 1 and lc_a != 0):
        return lc_a
    
    # If lc_a given as zero => lets FaceCenteredCubic get latticeconstant from ase
    elif(lc_a == 0):
        lc_a = None
        return lc_a
    
    # Returns a dict with all the values
    return lc_constants

""" Takes out the information from the ordered dictionary and creates an atomobject """
def from_dictionary_to_atoms(dictionary, symbol, Size_X, Size_Y, Size_Z):

    # Returns the chemical formula which is needed when creating the atoms object.
    chemical_formula_sum = str((dictionary[symbol])['_chemical_formula_sum'])
    
    # Sometimes the chemical formula have spaces in between the elements. This function removes the spaces. 
    chemical_formula = chemical_formula_sum.replace(" ","")
    
    # Lattice constants
    a = float((dictionary[symbol])['_cell_length_a'])
    b = float((dictionary[symbol])['_cell_length_b'])
    c = float((dictionary[symbol])['_cell_length_c'])
    
    # Angles for the unit cell
    alpha = float((dictionary[symbol])['_cell_angle_alpha'])
    beta = float((dictionary[symbol])['_cell_angle_beta'])
    gamma = float((dictionary[symbol])['_cell_angle_gamma'])
    
    # Nr of atoms in the structural chemical formula
    nr_of_atoms = len(((dictionary[symbol])['_atom_site_fract_x']))

    # Extracts the fractional coordinates that is then multipled with corresponding lattice constant
    # This will return the posistions between the atoms based on the structural chemical formula
    position = []
    for i in range(nr_of_atoms):
        x = float(((dictionary[symbol])['_atom_site_fract_x'])[i])*a
        y = float(((dictionary[symbol])['_atom_site_fract_y'])[i])*b
        z = float(((dictionary[symbol])['_atom_site_fract_z'])[i])*c
        pos = (x,y,z)
        position.append(pos)        

    # Creates the atomobject
    atoms =Atoms(symbols= chemical_formula,
                positions=position,
                cell=[a, b, c, alpha, beta, gamma])

    # Generates a supercell
    atoms = atoms*(Size_X,Size_Y,Size_Z)
    
    return atoms