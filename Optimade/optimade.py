import json
from ase.atoms import *

def translate_to_optimade(atomobj,MSD):
    Id = "blablabla" #Hitta på att få unika ID men samma för varje sample. Troligen en egen funktion
    #elements = 
    #nelements = 
    # element_ratios = 
    # Chemical_formula_anonymous =
    dimension_types = [1,1,1]
    nperiodic_dimensions = 3
    # lattice_vectrs =
    cartesian_site_positions = atomobj.get_positions()
    nsites = len(atomobj)
    species_at_sites = atomobj.get_chemical_symbols()
    # species =
    # structure_features = [] #oklart! 
    group1_Mean_Self_cofficient = MSD #Detta ska in via prefix kvantiter. Hur man gör detta har jag ingen aning om. 
    return None
