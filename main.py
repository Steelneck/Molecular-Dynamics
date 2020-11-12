"""Demonstrates molecular dynamics with constant energy."""

from init_values import *

def main():
    
    # Initiate the crystal based on the chosen variables
    # This will eventually become "Initiate the system" => system depends on user's choice
    init()
    
    # We want to run MD with constant energy using the VelocityVerlet algorithm.
    dyn = VelocityVerlet(atoms, 5 * units.fs)  # 5 fs time step.

    # Now run the dynamics
    dyn.attach(printenergy, interval=3)
    printenergy() # Use it for now to check results in terminal
    dyn.run(20)

if __name__ == "__main__":
    main()