from plot_functions import plot_prop_vs_time
from plot_functions import plot_prop_per_simulation
from plot_functions import plot_prop_all_elements
from plot_functions import hist_prop_per_simulation

"""
Script used to run plot functions. What plot will be shown depends on the input from plot_input.json

Available databases: "Materials Project" or "ASE"
Availablable integrators: "Velocity-Verlet" or "Langevin"
Available potentials: "EMT", "Leonnard-Jones" or  specific key for chosen kim potential as a string.
Examples:
"EAM_Dynamo_AcklandTichyVitek_1987_Cu__MO_179025990738_005"
"EAM_Dynamo_AcklandTichyVitek_1987_Ni__MO_977363131043_005"
"AM_Dynamo_AdamsFoilesWolfer_1989Universal6_Pd__MO_169076431435_000"
Available elements: Check sim_properties.csv columnt "Elements" to see which elements have been simulated.
For EMT: Cu, Ni, Au, Ag, Al, Pd, Pt (H, C, N, O only for fun appearently).
Available properties for plotting: "MSD", "D", "L", "SHC", "Int_T", "E_coh", "Int_P", "Bulk_mod".
Structures: Look in the sim_properties.csv too see which structures are saved.
"""

#Run plot functions. If multiple are active they will appear one after another.

plot_prop_vs_time("plot_input.json")
plot_prop_per_simulation("plot_input.json")
plot_prop_all_elements("plot_input.json")
hist_prop_per_simulation("plot_input.json")
