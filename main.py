"""Demonstrates molecular dynamics with constant energy."""

from Init.init_values import *

def main():
    
    # Initiate the crystal based on the chosen variables
    # This will eventually become "Initiate the system" => system depends on user's choice
    atoms = init()
    
    # We want to run MD with constant energy using the VelocityVerlet algorithm.
    dyn = VelocityVerlet(atoms, 5 * units.fs)  # 5 fs time step.
    
    """ Keep this for now to see changes in energy in terminal when changing lattice.
    This was a quick fix to solve the problem of attaching it to dyn """
    def printenergy(a=atoms):  # store a reference to atoms in the definition.
        """Function to print the potential, kinetic and total energy."""
        epot = atoms.get_potential_energy() / len(atoms)
        ekin = atoms.get_kinetic_energy() / len(atoms)
        print('Energy per atom: Epot = %.3feV  Ekin = %.3feV (T=%3.0fK)  '
        'Etot = %.3feV' % (epot, ekin, ekin / (1.5 * units.kB), epot + ekin))

    # Now run the dynamics
    dyn.attach(printenergy, interval=3)
    printenergy # Use it for now to check results in terminal
    dyn.run(20)

if __name__ == "__main__":
    main()
