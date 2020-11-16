"""Calculator for time averaged Mean Square Displacement"""
"""Unclear if MSD should be calculated from t after equilibrium or from equilibrium to t=0"""

from asap3 import Trajectory, FullNeighborList

import numpy as np
import tkinter

def MSD_calc(a, t, eq_list):
    try:
        a.get_masses() # Tries if the attribute exists, skips except if it does
    except AttributeError:
        print("You have not entered a valid system.") # Message for user if no attribute
        return False # Ends the function
    """Function to calculate and print the time average of the mean square displacement (MSD) at time t."""
    traj_MSD = Trajectory("atoms.traj")  #time evolution of trajectory
    time = len(traj_MSD)-t
    pos_eq = eq_list[-1].get_positions() #position of atoms when system has reached equilibrium
    pos_t = eq_list[t].get_positions() #position of atoms at time t
    diff = pos_t - pos_eq #displacement of all atoms from equilibrium to time t as a vector
    diff_sq = np.absolute(diff)**2 
    MSD = np.sum(diff_sq)/len(a) #Time averaged mean square displacement.
    print("MSD = ", MSD)
    return(MSD)

def Self_diffuse(MSD, t, eq_list):
    """Function that  calculates the self-diffusion coefficient (D) at time t, based on the value of the mean square displacement."""
    traj_MSD = Trajectory("atoms.traj")  
    time = len(traj_MSD)-t
    D = MSD/(6*time) #How to connect mean squre displacement to self-diffusion coefficient.
    open("atoms.traj", "w").close()
    print("D = ", D)
    return(D)

def Lindemann(a, MSD):
     nblist = FullNeighborList(3.5, a).get_neighbors(1, -1) #Returns 3 lists containing information about nearest neighbors. 3rd list is the square of the distance to the neighbors.
     d = np.sqrt(np.amin(nblist[2])) #distance to the nearest neighbor. Takes the minimum value of nblist.
     L = MSD/d #Lindemann criterion. Expect melting when L>0.1
     if L > 0.1:
         #print("MELTING!")
         print(L)
         return True
     else:
         #print("Not melting.")
         print(L)
         return False




