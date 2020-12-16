import os, csv, json, operator

def write_json_to_csv(path):
    files = os.listdir(path)
    for file in files:
        if file.endswith(".json"): #Only want to run through .json-files
            with open(path + file, "r") as json_file:
                data = json.load(json_file)
                for sim_in in data["Simulation input"]:
                    Chem_form = sim_in["Chemical formula"]
                    Database = sim_in["Database"]
                    Potential = sim_in["Potential"]
                    Integrator = sim_in["Integrator"]
                    Struct = sim_in["Lattice structure"]
                    Init_T = sim_in["Initial Temperature"]
                    Run_time = sim_in["Run time"]
                    Volume = sim_in["Volume"]
                for sim_out in data["Simulation output"]:
                    MSD = sim_out["Mean Square Displacement"]
                    D = sim_out["Self diffusion coefficient"]
                    L = sim_out["Lindemann Criterion"]
                    SHC = sim_out["Specific Heat"]
                    Int_T = sim_out["Internal Temperature"]
                    E_coh = sim_out["Cohesive Energy"]
                    Int_P = sim_out["Internal Pressure"]
                    Bulk_mod = sim_out["Bulk_modulus"]
                    Lattice_constant = sim_out["Lattice constant"]

            csvfilepath = "Visualization/sim_properties.csv"
            with open(csvfilepath, "a", newline = "") as csvfile:
                fieldnames = ["Element",
                              "Database",
                              "Potential",
                              "Integrator",
                              "Struct",
                              "Init_T",
                              "Run_time",
                              "Vol",
                              "MSD",
                              "D",
                              "L",
                              "SHC",
                              "Int_T",
                              "E_coh",
                              "Int_P",
                              "Bulk_mod",
                              "Lattice_constant"]
                
                writer = csv.DictWriter(csvfile, fieldnames = fieldnames, delimiter=";")
                if os.path.getsize(csvfilepath) == 0: #Write header if file is empty
                    writer.writeheader()
                writer.writerow({"Element" : Chem_form,
                                 "Database" : Database,
                                 "Potential" : Potential,
                                 "Integrator" : Integrator,
                                 "Struct" : Struct,
                                 "Init_T" : Init_T,
                                 "Run_time" : Run_time,
                                 "Vol" : Volume,
                                 "MSD" : MSD,
                                 "D" : D,
                                 "L" : L,
                                 "SHC" : SHC,
                                 "Int_T" : Int_T,
                                 "E_coh" : E_coh,
                                 "Int_P" : Int_P,
                                 "Bulk_mod" : Bulk_mod,
                                 "Lattice_constant": Lattice_constant})

write_json_to_csv("Json/")
