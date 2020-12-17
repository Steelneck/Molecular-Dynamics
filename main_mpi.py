"""Demonstrates molecular dynamics with constant energy."""

from Init.init_mpi_simulation import *
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
            Calculator_Check = row["Calculator"]
            Algorithm_Check = row["Algorithm"]
            KIM_potential = row["KIM_potential"]
            Symbol = row["Symbol"]
            Defect_Check = row["Defect"]
            Impurity_ele_list=row["Impurity_ele_list"]
            Temperature = row["Temperature"]
            Steps = row["Steps"]
            Interval = row["Interval"]
            Size_X = row["Size_X"]
            Size_Y = row["Size_Y"]
            Size_Z = row["Size_Z"]
            API_Key = row["API_Key"]
            PBC = row["PBC"]
            Directions = row["Directions"]
            Miller = row["Miller"]
            lc_a=row["lc_a"]
            lc_b=row["lc_b"]
            lc_c=row["lc_c"]
            lc_alpha=row["lc_alpha"]
            lc_beta=row["lc_beta"]
            lc_gamma=row["lc_gamma"]

        simulation_mpi(Calculator_Check, 
                        Algorithm_Check, 
                        KIM_potential,
                        Symbol,
                        Defect_Check,
                        Impurity_ele_list,
                        Temperature,
                        Steps,
                        Interval,
                        Size_X,
                        Size_Y,
                        Size_Z,
                        API_Key,
                        PBC,
                        Directions,
                        Miller,
                        lc_a,
                        lc_b,
                        lc_c,
                        lc_alpha,
                        lc_beta,
                        lc_gamma)

if __name__ == "__main__":
    main()