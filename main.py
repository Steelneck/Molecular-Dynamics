"""Demonstrates molecular dynamics with constant energy."""

from Init.init_simulation import *

def main():  
    
    simulation(EMT(),'Cu')
    simulation(EMT(), "Ni")
    simulation(EMT(), "Pd")
    
if __name__ == "__main__":
    main()
