"""Demonstrates molecular dynamics with constant energy."""

from Init.init_simulation import *

def main():  
    
<<<<<<< HEAD
    simulation()
    
=======
    for atomobj in atoms:
        
        # Describe the interatomic interactions with the Effective Medium Theory
        # (Note: Create a higher ordet function)
        atomobj.calc = Calculator
        
        # We want to run MD with constant energy using the VelocityVerlet algorithm.
        dyn = VelocityVerlet(atomobj, 5*units.fs)  # 5 fs time step.


        traj = Trajectory("atoms.traj", "w", atomobj)

        dyn.attach(traj.write, interval)
        dyn.run(steps)
        
        traj = Trajectory("atoms.traj")
        traj_eq = Trajectory("atoms_eq.traj", "w", atomobj)
        
        calc.eq_traj(atomobj, traj, traj_eq, Size_X * Size_Y * Size_Z)#Creates new .traj-file containing trajectory post equilibrium.
        if os.path.getsize("atoms_eq.traj") != 0: #If-statement that checks if we ever reached equilibrium. Returns a message if the traj-file is empty, otherwise does calculations.
            traj_eq = Trajectory("atoms_eq.traj")
            #If-statement that checks if we ever reached equilibrium.
            MSD = calc.MSD_calc(atomobj, traj_eq, timeStepIndex)
            D = calc.Self_diffuse(traj_eq, MSD)
            L = calc.Lindemann(traj_eq, MSD)
            SHC = calc.Specific_Heat(atomobj, traj_eq)

            # Internal temperature of the system
            internalTemperature = calc.internal_temperature(atomobj, traj_eq, timeStepIndex)

        else:
            print("System never reached equilibrium. No calculations are possible.")

>>>>>>> feature_Materials_project
if __name__ == "__main__":
    main()
