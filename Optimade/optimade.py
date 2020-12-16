import json
import hashlib
from ase.atoms import *
from pymatgen import Composition
from ase import cell
from datetime import datetime, timezone
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

def translate_to_optimade(atomobj, meansSquareDisplacement, selfDiffusionCoffecient, lindemann , specificHeatCapacity, 
                                    internalTemperature, cohesiveEnergy, internalPressure, bulkModulus):
    #Creates a cell from the atomobject
    cell = atomobj.cell

    nsites = len(atomobj)
    # Extract all data from atoms object
    Comp = Composition(atomobj.get_chemical_formula())
    id = 0                                                  # Init an id, changes below.

    # Creates a timestamp in iso format 
    local_time = datetime.now(timezone.utc).astimezone()
    local_time.isoformat()

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
    
    nelements = len(elements)                        # Number of _different_ elements
    Pretty_formula = Comp.reduced_formula
    Chemical_formula_anonymous = Comp.anonymized_formula
    Chemical_formula_descriptive = atomobj.get_chemical_formula()
    Chemical_formula_reduced = Pretty_formula
    dimension_types = [1,1,1]
    nperiodic_dimensions = 3
    lattice_vectors = cell[:]
    cartesian_site_positions = atomobj.get_positions() 
    species_at_sites = atomobj.get_chemical_symbols()
    # species = 
     
    
    # Init a dictionary
    data_dict = {}
    # --- Add data with all required keys
    data_dict["_id"] = id
    data_dict["mean_square_displacement"] = meansSquareDisplacement
    data_dict["self_diffusion"] = selfDiffusionCoffecient
    data_dict["lindemann"] = lindemann
    data_dict["specific_heat"] = specificHeatCapacity
    data_dict["internal_temperature"] = internalTemperature
    data_dict["cohesive_energy"] = cohesiveEnergy
    data_dict["internal_pressure"] = internalPressure
    data_dict["bul_modulus"] = bulkModulus
    data_dict["cartesian_site_positions"] = cartesian_site_positions.tolist()
    data_dict["dimension_types"] = dimension_types
    data_dict["nperiodic_dimensions"] = nperiodic_dimensions
    data_dict["elements"] = elements
    data_dict["element_ratios"] = element_ratios
    data_dict["formula_anonymous"] = Chemical_formula_anonymous
    data_dict["last_modified"] = { "$date" : local_time.isoformat() + "Z"}
    data_dict["lattice_vectors"] = lattice_vectors.tolist()
    data_dict["nelements"] = nelements
    data_dict["nsites"] = nsites
    data_dict["chemical_formula_descriptive"] = Chemical_formula_descriptive
    data_dict["chemical_formula_reduced"] = Chemical_formula_reduced
    data_dict["structure_feature"] : []
    #print(element_ratios)


    # Make JSON object
    data_json = json.dumps(data_dict)
        # Create unique id from contents
    unique_id = hashlib.md5(data_json.encode("utf-8")).hexdigest()
    data_dict["_id"] = {"oid" : unique_id}
    # Add id to data object and recreate a json object to write to destinationfile
    with open("Optimade/data_format_optimade.json", 'w') as json_data_file:
        json.dump(data_dict, json_data_file, indent=2)

    print("Data for optimade is temporarly stored in 'Optimade/data_format_optimade.json'.")


"""
[
    {
      "_id": {
        "$oid": "5cfb441f053b174410700d02"
      },
      "mean_square_displacement": 12,
      "self_diffusion": 12,
      "lindemann": "meling",
      "specific_heat": 12,
      "internal_temperature": 12,
      "cohesive_energy": 12,
      "internal_pressure": 12,
      "bul_modulus": 12,
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
          1.666666666666666,
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
          0.28
        ]
      ],
      "nelements": 1,
      "nsites": 1,
      "chemical_formula_descriptive" : "CuCu",
      "chemical_formula_reduced" : "CuH",
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
  }
  ]
"""