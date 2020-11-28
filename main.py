"""Demonstrates molecular dynamics with constant energy."""

<<<<<<< HEAD
from Init.init_values import *
from tkinter import *
=======
from Init.init2 import *
>>>>>>> Created a new init file for materials project. Can now extract information and create atoms objects. Need to do more research on this
import os
from tkinter import *
import Calculations.calculations as calc
from asap3 import Trajectory

def main():  
    # Initiate the crystal based on the chosen variables
    # This will eventually become "Initiate the system" => system depends on user's choice
<<<<<<< HEAD
    atoms = init()
=======
    
    atoms = init_MP()
>>>>>>> Created a new init file for materials project. Can now extract information and create atoms objects. Need to do more research on this
    
    
    # We want to run MD with constant energy using the VelocityVerlet algorithm.
    dyn = VelocityVerlet(atoms, 5*units.fs)  # 5 fs time step.

<<<<<<< HEAD
    traj = Trajectory("atoms.traj", "w", atoms)

    dyn.attach(traj.write, interval=10)
    dyn.run(steps=200)
    
    traj = Trajectory("atoms.traj")
    traj_eq = Trajectory("atoms_eq.traj", "w", atoms)
    
    calc.eq_traj(atoms, traj, traj_eq, Size_X * Size_Y * Size_Z)#Creates new .traj-file containing trajectory post equilibrium.
    if os.path.getsize("atoms_eq.traj") != 0: #If-statement that checks if we ever reached equilibrium. Returns a message if the traj-file is empty, otherwise does calculations.
        traj_eq = Trajectory("atoms_eq.traj")
        #If-statement that checks if we ever reached equilibrium.
        MSD = calc.MSD_calc(atoms, 10)
        D = calc.Self_diffuse(MSD, 10)
        L = calc.Lindemann(atoms, MSD)
        SHC = calc.Specific_Heat(atoms)

        # Internal temperature of the system
        internalTemperature = calc.internal_temperature(atoms, traj_eq, 10)
=======
        traj = Trajectory("atoms.traj", "w", atomobj)
        
        dyn.attach(traj.write, interval=10)
        dyn.run(100)
>>>>>>> Created a new init file for materials project. Can now extract information and create atoms objects. Need to do more research on this
        
        # Debye temperature of the system
        #debyeTemperature = calc.debye_temperature(atoms, "atoms_eq.traj", 10, EMT())

        internalPressure = calc.calc_internal_pressure(atoms, traj_eq, Size_X * Size_Y * Size_Z)
    else:
        print("System never reached equilibrium. No calculations are possible.")

if __name__ == "__main__":
    main()
