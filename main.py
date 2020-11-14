from ase.lattice.cubic import FaceCenteredCubic
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
                          size=(1, 1, 1),
                          pbc=True)
    

    # Describe the interatomic interactions with the Effective Medium Theory
    calc = KIM("LJ_ElliottAkerson_2015_Universal__MO_959249795837_003")
    atoms.calc = calc

    # Set the momenta corresponding to T=300K
    MaxwellBoltzmannDistribution(atoms, 300 * units.kB)

    # We want to run MD with constant energy using the VelocityVerlet algorithm.
    dyn = VelocityVerlet(atoms, 5 * units.fs)  # 5 fs time step.
    # dyn = Langevin(atoms, 5 * units.fs, units.kB * 300, 0.002)

    traj = Trajectory("atoms.traj", "w", atoms)
    dyn.attach(traj.write, interval=1)
    dyn.run(20)
    
    MSD = MSD_calc(atoms, 10)
    D = Self_diffuse(MSD, 10)
    L = Lindemann(atoms, MSD)

if __name__ == "__main__":
    main()
