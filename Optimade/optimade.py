import json
import hashlib
from ase.atoms import *
from pymatgen import Composition
from ase import cell
from datetime import datetime
import glob, os
import shutil
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
    # Create unique id
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
      raise ValueError("Could not concatenate data files for Optimade. \n Please provide the run_id as a string. Provided type where:", type(run_id))
  except ValueError as e:
    print(e)
    print("Saving concatenated file with run_id: TEMPORARY_ID. If you don't change filename until next run, this file will be overwritten.")
    run_id = "TEMPORARY_ID"
    pass
  
  # Optimade directory
  optimadeDir = "Optimade/"

  # Filter out all .json files in Optimade directory
  allDataFiles = []
  for file in os.listdir(optimadeDir):
    if file.endswith(".json"):
      allDataFiles.append(file)

  if allDataFiles:                                                             # Python style to check non-empty list
    # Extract existing concatenated files and remove them from allDataFiles list. We don't want to append those.
    concatenatedFiles = [conc for conc in allDataFiles if "concatenated_data_for_optimade_" in conc]
    for concFile in concatenatedFiles:
      allDataFiles.remove(concFile)

    # --- Create targetfile with concatenated data
    fileFormat = ".json"
    targetFileName = "concatenated_data_for_optimade_" + str(run_id)

    # Generate file extension (a number) if file already exists.
    fileExtension = ""
    if os.path.exists(optimadeDir + targetFileName + fileFormat):
      print("Concatenated file with run id:", str(run_id), "already exists. Appending file extension...")
      fileExtensionNumber = 1
      # Loop until we find an unused file extension
      while os.path.exists(optimadeDir + targetFileName + fileExtension + fileFormat):
        fileExtension = "_" + str(fileExtensionNumber)
        fileExtensionNumber += 1
    
    targetFile = optimadeDir + targetFileName + fileExtension + fileFormat      # Set targetFileJSON with generated extension
    # ---

    # Create a folder to store all datafiles if folder does not exist already
    sourceDataFolder = optimadeDir + "Source_data_" + str(run_id) + fileExtension
    if not os.path.exists(sourceDataFolder):
      os.mkdir(sourceDataFolder)

    allDataList = []                                                            # Init a final list that will be written to concatenated file
    for file in allDataFiles:                                                   # Open all files and append json data to final list.
      with open(optimadeDir + file) as json_data_file:
        json_data = json.load(json_data_file)         
        allDataList.append(json_data[0])                                        # Dict is stored in a list in each json file – extract first (and only) item
        shutil.move(optimadeDir + file, sourceDataFolder + "/" + file)          # Move file to subdirectory to be stored.
      json_data_file.close()

    print("Concatenated file is stored in: " + targetFile)
    with open(targetFile, 'w') as targetFileJSON:
      json.dump(allDataList, targetFileJSON, indent=2)                          # Finally write concatenated data to target file
    targetFileJSON.close()
  else: 
    print("No JSON data files in Optimade folder.")

    