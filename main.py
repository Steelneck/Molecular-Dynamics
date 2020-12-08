"""Demonstrates molecular dynamics with constant energy."""

from Init.init_simulation import *
import json

def main():  

    with open('User_Input.json') as json_file:
        Input = json.load(json_file)
        for row in Input['Data']:
            Name = row["Name"]
            EMT = row["EMT"]
            openKIM = row["openKIM"]
            KIM_potential = row["KIM_potential"]
            ASE = row["ASE"]
            Materials_project = row["Materials_project"]
            Symbol = row["Symbol"]
            Temperature = row["Temperature"]
            Steps = row["Steps"]
            Interval = row["Interval"]
            Size_X = row["Size_X"]
            Size_Y = row["Size_Y"]
            Size_Z = row["Size_Z"]
            simulation(Name, EMT,openKIM,KIM_potential,ASE, Materials_project,
                        Symbol, Temperature, Steps, Interval,
                        Size_X, Size_Y, Size_Z)

    
if __name__ == "__main__":
    main()
