import json
import glob, os
import shutil

def main():
    print("Starting data transfer to optimade repository. Make sure that all desired concatenated data files are in directory '/Optimade'. \n")

    # Ask user to provide full path to target json file, that will be uploaded to optimade.
    optimadePath = input("Please enter absolute path to optimade-python-tools repository server structure data file. Example: <my_path>/optimade-python-tools/optimade/server/data/<my_structures_file>.json \n")
    # Remove possebly accidental spaces
    optimadePath = optimadePath.strip(" ")                                          
    while not os.path.exists(optimadePath) or ".json" not in optimadePath:
        print("Provided optimade path:", optimadePath, "does not exist or is not a .json file. Please try again.")
        optimadePath = input("Please enter absolute path to optimade-python-tools repository server structure data file. Example: <my_path>/optimade-python-tools/optimade/server/data/<my_structures_file>.json \n")


    structureData = []
    try:
        with open(optimadePath) as structureDataFile:
            structureData = json.load(structureDataFile)
        structureDataFile.close()
    except IOError as e:
        print(e)
        return None
    
    # Ask user to confirm
    response = input("\nAll data in the files starting with 'concatenated_data_for_optimade_' in directory '/Optimade' will be transfered to provided JSON file: '" + optimadePath + "'. Make sure only desired runs are in /Optimade folder and target location is correct. Approve/deny transfer with yes/no:\n")
    # Remove possebly accidental spaces
    response = response.strip(" ")
    
    # Check that response is either yes, no, y, n. Could not come up with a way not to copy code as following, without creating loop holes. 
    responseCheck = True
    while responseCheck:
        cwd = os.getcwd()                                                               # Current working directory to change back to below
        if response == "yes" or response == "y":
            os.chdir(os.path.abspath("Optimade"))                                           # Changes directory to Optimade, holds for rest of function
            allDataFiles = glob.glob("*.json")                                              # Extract all .json files in Optimade folder
            
            # Create a folder to store all transfered datafiles if folder does not exist already
            transferHistoryFolder = "Transfer_History"
            if not os.path.exists(transferHistoryFolder):
                os.mkdir(transferHistoryFolder)

            # Extract existing concatenated files and remove them from allDataFiles list. We don't want to append those.
            concatenatedFiles = [conc for conc in allDataFiles if "concatenated_data_for_optimade_" in conc]
            print("\nThe following data files will be merged into: '" + optimadePath + "'. You can find the transfer history in folder '" + transferHistoryFolder + "'.")
            for concFile in concatenatedFiles:
                print("   ", concFile)
                # Open each file and pack into target JSON file
                with open(concFile) as sourceDataFile:
                    sourceData = json.load(sourceDataFile)
                
                shutil.move(concFile, transferHistoryFolder + "/" + concFile)               # Move file to subdirectory Transfer_History to be stored as reference.
                for sourceDict in sourceData:
                    structureData.append(sourceDict)                                         # Structure in JSON files are list(dict(), ...) hence first element is our source dict
                sourceDataFile.close()
            

            # Open target file in write mode and overwrite with updated data
            with open(optimadePath, 'w') as structureDataFile:
                json.dump(structureData, structureDataFile, indent=2)                       # structureData should be on list format still.
            
            os.chdir(cwd)                                                                   # Change back to original directory
            print("\nTransfer complete, you may now run local optimade server to check success.")
            break

        elif response == "no" or response == "n":
            os.chdir(cwd)                                                                   # Change back to original directory
            print("\nStopping transfer, re-run program if you want to try again.")
            break
        else:
            print("Please type yes/no.")
            response = input("All data in the files starting with 'concatenated_data_for_optimade_' in directory '/Optimade' will be transfered to provided JSON file: '" + optimadePath + "'. Make sure only desired runs are in /Optimade folder and target location is correct. Approve/deny transfer with yes/no:\n")
            continue


if __name__ == "__main__":
    main()