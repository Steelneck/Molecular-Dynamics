from ase.lattice.cubic import FaceCenteredCubic
from ase.lattice.cubic import BodyCenteredCubic
from ase.lattice.cubic import SimpleCubic
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.md.verlet import VelocityVerlet
from ase.lattice.bravais import Lattice
from ase.atoms import Atoms
from ase.md.langevin import Langevin
from ase import units
from ase.calculators.kim.kim import KIM

from asap3 import EMT, Trajectory
from asap3.io.trajectory import *

from Calculations.MSD import MSD_calc
from Calculations.MSD import Lindemann
from Calculations.MSD import Self_diffuse


def main():
    atoms = FaceCenteredCubic(directions=[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
                          symbol="Cu",
                          size=(5, 5, 5),
                          pbc=True)
    

    # Describe the interatomic interactions with the Effective Medium Theory
    #calc = KIM("LJ_ElliottAkerson_2015_Universal__MO_959249795837_003")
    #atoms.calc = calc
    atoms.calc = EMT()

    # Set the momenta corresponding to T=300K
    MaxwellBoltzmannDistribution(atoms, 300 * units.kB)

    # We want to run MD with constant energy using the VelocityVerlet algorithm.
    dyn = VelocityVerlet(atoms, 5 * units.fs)  # 5 fs time step.
    # dyn = Langevin(atoms, 5 * units.fs, units.kB * 300, 0.002)
    
    traj = Trajectory("atoms.traj", "w", atoms)
    
    dyn.attach(traj.write, interval=10)
    dyn.run(200)

    traj_MSD = Trajectory("atoms.traj")
    n = 1
    atoms_eq = []
    print(traj_MSD)
    while n < len(traj_MSD):
        ediff = abs((traj_MSD[n].get_potential_energy() / len(atoms)) -
               (traj_MSD[n].get_kinetic_energy() / len(atoms))) #epot-ekin
        if ediff < 0.005:
            atoms_eq.append(traj_MSD[n])
        n += 1
    print(len(atoms_eq))
    MSD = MSD_calc(atoms, 10, atoms_eq)
    D = Self_diffuse(MSD, 10, atoms_eq)
    L = Lindemann(atoms, MSD)

    open("atoms.traj", "w").close()
    
if __name__ == "__main__":
    main()
