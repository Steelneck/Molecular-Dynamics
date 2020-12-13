"""Demonstrates molecular dynamics with constant energy."""

from Init.init_super_simulation import *

def main():  
    
    super_simulation(EMT(),'Cu')
    super_simulation(EMT(),'Pd')
    super_simulation(EMT(),'Ni')
    super_simulation(EMT(),'Ag')
    
if __name__ == "__main__":
    main()
