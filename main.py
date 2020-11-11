from ase.lattice.cubic import FaceCenteredCubic
from ase import Atoms
from Calculations.calculations import Specific_Heat

def main():
    fcclattice = FaceCenteredCubic(directions=[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
                              symbol='Cu',
                              size=(1, 1, 1),
                              pbc=True)
    print(Specific_Heat(fcclattice))

    atoms = Atoms('CO', positions=[(0, 0, 0),(0,0,2)])
    print(Specific_Heat(atoms))

if __name__ == "__main__":
    main()
    