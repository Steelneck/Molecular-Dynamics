import json
import hashlib
from ase.atoms import *
from pymatgen import Composition
from ase import cell
from datetime import datetime
import glob, os
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
    local_time = datetime.now()
    local_time.isoformat()

    """
    One method if formula is depricated as stated in ase documentation, however the method below is cleaner.
    Also makes 'elements' sync with 'element_ratios'. 
    elements = []
    for atomSymbol in atomobj.get_chemical_symbols():
        if atomSymbol not in elements:
            elements.append(atomSymbol)
    """
    symbols_count_item_list = sorted(atomobj.symbols.formula.count().items())   # Sorted list of items = [key, value]
    elements = []
    element_ratios = []
    for item in symbols_count_item_list:                # Key at item[0], value at item[1]
      elements.append(item[0])
      element_ratios.append(float(item[1]/nsites))      # Ratio of each element in structure. Float-cast to ensure floats
    
    nelements = len(elements)                           # Number of _different_ elements
    Chemical_formula_anonymous = Comp.anonymized_formula
    Chemical_formula_descriptive = atomobj.get_chemical_formula()
    Chemical_formula_reduced = atomobj.get_chemical_formula()
    dimension_types = [1,1,1]
    nperiodic_dimensions = 3
    element_amount_dict = Comp.get_el_amt_dict()
    lattice_vectors = cell[:]
    cartesian_site_positions = atomobj.get_positions() 
    species_at_sites = atomobj.get_chemical_symbols()
    species_list = []
    
    for x in elements:
      species_list.append({"chemical_symbols" : [x], "concentration" : [element_amount_dict[x]], "name": x})

    # Sets the disorder flag if there are more than 1 element present
    
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
    data_dict["bulk_modulus"] = bulkModulus
    data_dict["cartesian_site_positions"] = cartesian_site_positions.tolist()
    data_dict["dimension_types"] = dimension_types
    data_dict["nperiodic_dimensions"] = nperiodic_dimensions
    data_dict["elements"] = elements
    data_dict["elements_ratios"] = element_ratios
    data_dict["formula_anonymous"] = Chemical_formula_anonymous
    data_dict["last_modified"] = { "$date" : local_time.isoformat() + "Z"}
    data_dict["lattice_vectors"] = lattice_vectors.tolist()
    data_dict["nelements"] = nelements
    data_dict["nsites"] = nsites
    data_dict["chemical_formula_descriptive"] = Chemical_formula_descriptive
    data_dict["chemical_formula_reduced"] = Chemical_formula_reduced
    data_dict["species"] = species_list
    data_dict["species_at_sites"] = species_at_sites
    data_dict["structure_features"] = []
    data_dict["task_id"] = id                                                                       # Should these have unique id also?
    data_dict["relationships"] = {"references" : {"data" : [{"type" : "references", "id" : id}]}}   # Should these have unique id also?


    # Make JSON object
    data_json = json.dumps(data_dict)
        # Create unique id from contents
    unique_id = hashlib.md5(data_json.encode("utf-8")).hexdigest()
    data_dict["_id"] = unique_id
    # Add id to data object and recreate a json object to write to destinationfile
    fileName = "Optimade/data_for_optimade_" + str(unique_id) + ".json"
    with open(fileName, 'w') as json_data_file:
        json.dump([data_dict], json_data_file, indent=2)

    json_data_file.close()

    print("Data for optimade is temporarly stored in " + fileName + ".")


def concatenateOptimadeDataFiles(run_id):
  try:
    if type(run_id) != str:
      raise ValueError("Please provide the run_id as a string. Provided type where:", type(run_id))
  except ValueError as e:
    print(e)
    pass

  os.chdir(os.path.abspath("Optimade"))             # Changes directory to Optimade, holds for rest of function
  allDataFiles = glob.glob("*.json")                # Extract all .json files in Optimade folder

  # Extract existing concatenated files and remove them from allDataFiles list. We don't want to append those.
  concatenatedFiles = [conc for conc in allDataFiles if "concatenated_data_for_optimade_" in conc]
  for concFile in concatenatedFiles:
    allDataFiles.remove(concFile)

  allDataList = []                                  # Init a final list that will be written to concatenated file

  for file in allDataFiles:                         # Open all files and append json data to final list.
    with open(file) as json_data_file:
      json_data = json.load(json_data_file)         
      allDataList.append(json_data[0])              # Dict is stored in a list in each json file – extract first (and only) item
    json_data_file.close()

  targetFileName = "concatenated_data_for_optimade_" + str(run_id) + ".json"
  with open(targetFileName, 'w') as targetFile:
    json.dump(allDataList, targetFile, indent=2)    # Finally write concatenated data to target file
  targetFile.close()

  # Might want to move all files that where concatenated to a subfolder here. 

    

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
      "bulk_modulus": 12,
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
              "id": "Group1",

            }
          
          ]
      
      }
    }
  }
  ]
"""