import json
import hashlib
from ase.atoms import *
from ase.formula import Formula

"""
This is the interface to convert the simulation output to optimade format.

Required root keys are: 
_id
type
last_modified
species
    name
    chemical_symbols
    concentration


"""

def translate_to_optimade(atomobj,MSD):

    # Extract all data from atoms object
    id = 0                                                  # Init an id, changes below.
    nsites = len(atomobj)
    """
    One method if formula is depricated as stated in ase documentation, however the method below is cleaner.
    Also makes 'elements' sync with 'element_ratios'. 
    elements = []
    for atomSymbol in atomobj.get_chemical_symbols():
        if atomSymbol not in elements:
            elements.append(atomSymbol)
    """
    symbols_count_dict = atomobj.symbols.formula.count()
    elements = list(symbols_count_dict.keys())          # All elements as Chemical symbols in a list (only type)
    element_ratios = []
    for amount in symbols_count_dict.values():
        element_ratios.append(float(amount/nsites))     # Ratio of each element in structure. Float-cast to ensure floats
    
    ###!!! Verkar inte vara en grej i dokumentationen men den är i exemplet....
    # chemsys = "".join(elements)                         # All elements in a string
    nelements = len(elements)                           # Number of _different_ elements
   
    # Chemical_formula_anonymous =
    dimension_types = [1,1,1]
    nperiodic_dimensions = 3
    # lattice_vectrs =
    cartesian_site_positions = atomobj.get_positions()
    species_at_sites = atomobj.get_chemical_symbols()
    # species =
    # structure_features = [] #oklart! 
    group1_Mean_Self_cofficient = MSD #Detta ska in via prefix kvantiter. Hur man gör detta har jag ingen aning om. 
     
    # Init a dictionary
    data_dict = {}
    # --- Add data with all required keys
    data_dict["_id"] = id
    data_dict["elements"] = elements
    data_dict["nelements"] = nelements
    data_dict["nsites"] = nsites
    data_dict["elements"] = elements
    data_dict["element_ratios"] = element_ratios
    
    #print(element_ratios)

    # Make JSON object
    data_json = json.dumps(data_dict)

    # Create unique id from contents
    unique_id = hashlib.md5(data_json.encode("utf-8")).hexdigest()
    
    # Add id to data object and recreate a json object to write to destinationfile
    data_dict["_id"] = {"oid" : unique_id}
    with open("Optimade/data_format_optimade.json", 'w') as json_data_file:
        json.dump(data_dict, json_data_file, indent=2)

    print("Data for optimade is temporarly stored in 'Optimade/data_format_optimade.json'.")


"""
{
    "_id": {
      "$oid": "5cfb441f053b174410700d02"
    },
    "chemsys": "Ac",
    "cartesian_site_positions": [
      [
        0.17570227444196573,
        0.17570227444196573,
        0.17570227444196573
      ]
    ],
    "dimension_types": [
      1,
      1,
      1
    ],
    "nperiodic_dimensions": 3,
    "elements": [
      "Ac"
    ],
    "elements_ratios": [
      1.0
    ],
    "formula_anonymous": "A",
    "last_modified": {
      "$date": "2019-06-08T05:13:37.331Z"
    },
    "lattice_vectors": [
      [
        1.2503264826932692,
        0,
        0
      ],
      [
        0,
        9.888509716321765,
        0
      ],
      [
        0,
        0,
        0.2972637673241818
      ]
    ],
    "nelements": 1,
    "nsites": 1,
    "pretty_formula": "Ac",
    "species": [
      {
        "chemical_symbols": [
          "Ac"
        ],
        "concentration": [
          1.0
        ],
        "name": "Ac"
      }
    ],
    "species_at_sites": [
      "Ac"
    ],
    "structure_features": [],
    "task_id": "mpf_1",
    "relationships": {
      "references": {
        "data": [
          {
            "type": "references",
            "id": "dijkstra1968"
          }
        ]
      }
    }
"""