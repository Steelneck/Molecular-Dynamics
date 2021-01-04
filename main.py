"""Demonstrates molecular dynamics with constant energy."""

from Init.init_simulation import *
import json

""" Ignore the parallel rank deprecation warning """

def warn(*args, **kwargs):
    pass
import warnings
warnings.warn = warn

""" Main """

def main():  
    """ Takes all the data from User_Input.json """
    with open('User_Input.json') as json_file:
        Input = json.load(json_file)
        for row in Input['Data']:
            EMT_Check = row["EMT"]
            openKIM_Check = row["openKIM"]
            Lennard_Jones_Check = row["Lennard_Jones"]
            LJ_epsilon = row["LJ_epsilon"]
            LJ_sigma = row["LJ_sigma"]
            LJ_cutoff = row["LJ_cutoff"]
            Velocity_Verlet_Check = row["Velocity_Verlet"]
            Langevin_Check = row["Langevin"]
            Langevin_friction = row["Langevin_friction"]
            Time_step = row["Time_step"]
            KIM_potential = row["KIM_potential"]
            ASE = row["ASE"]
            Symbol = row["Symbol"]
            Materials_project = row["Materials_project"]
            API_Key = row["API_Key"]
            Criteria_list=row["Criteria_list"]
            Vacancy = row["Vacancy"]
            Impurity = row["Impurity"]
            Impurity_ele_list=row["Impurity_ele_list"]
            Temperature = row["Temperature"]
            Steps = row["Steps"]
            Interval = row["Interval"]
            Size_X = row["Size_X"]
            Size_Y = row["Size_Y"]
            Size_Z = row["Size_Z"]
            PBC = row["PBC"]
            Bravais_lattice = row["Bravais_lattice"]
            Directions = row["Directions"]
            Miller = row["Miller"]
            lc_a=row["lc_a"]
            lc_b=row["lc_b"]
            lc_c=row["lc_c"]
            lc_alpha=row["lc_alpha"]
            lc_beta=row["lc_beta"]
            lc_gamma=row["lc_gamma"]
            Run_Optimade=row["Run_Optimade"]
            Optimade_name = row["Optimade_name"]

            simulation(EMT_Check,openKIM_Check, Lennard_Jones_Check, LJ_epsilon,
                        LJ_sigma, LJ_cutoff, Velocity_Verlet_Check, 
                        Langevin_Check, Langevin_friction, Time_step, KIM_potential,
                        ASE, Symbol, Materials_project,API_Key,Criteria_list, 
                        Vacancy, Impurity, Impurity_ele_list,
                        Temperature, Steps, Interval,Size_X, Size_Y, Size_Z,
                        PBC, Bravais_lattice,Directions,Miller,
                        lc_a,lc_b,lc_c,lc_alpha,lc_beta,lc_gamma,Run_Optimade,Optimade_name)

if __name__ == "__main__":
    main()
