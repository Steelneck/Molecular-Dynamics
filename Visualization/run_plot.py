import csv, operator
import matplotlib.pyplot as plt
import numpy as np

def plot_prop_vs_time(csvFileName, prop_str):

    """
    Function which plots the time evolution of property specified by input prop_str. csvFileName is the .csv-file where the time-evolution is saved."
    """
    file = open(csvFileName, "r")
    reader = csv.DictReader(file, delimiter=";")
    count = 1
    x_list = []
    y_list = []
    #Go through all rows in the csvfile and write append to lists.
    for row in reader: 
        x_list.append(float(row["Time"])) #Uses data from column labeled "Time"
        y_list.append(float(row[prop_str])) #Uses data from column labeled prop_str.
        count += 1
        
    plt.scatter(x_list,y_list)

    #Adds trendline to the plot.
    z = np.polyfit(x_list, y_list, 1)
    p = np.poly1d(z)
    plt.plot(x_list, p(x_list), "r--")

    plt.title(prop_str + " vs time scatter plot")
    plt.xlabel("Time [fs]")
    plt.ylabel(prop_str)

    #Set interval for the y-axis to not have the plot zoomed in too much.
    plt.ylim(0,0.1)

    plt.show()

def plot_prop_per_simulation(atomName,
                             struct_str,
                             prop1_str,
                             prop2_str,
                             database,
                             integrator,
                             potential):
    
    file = open("sim_properties.csv", "r")
    reader = csv.DictReader(file, delimiter=";")
    count = 1
    x_list = []
    y_list = []
    for row in reader:
        if row["Database"] == database and row["Integrator"] == integrator:
            if row["Potential"] == potential:
                if row["Element"] == atomName and row["Struct"].startswith(struct_str):
                    x_list.append(float(row[prop1_str]))
                    y_list.append(float(row[prop2_str]))
                    count += 1             
    plt.scatter(x_list, y_list)

    #These lines of code are only needed if a trendline is requested.
    #z = np.polyfit(x_list, y_list, 1)
    #p = np.poly1d(z)
    #plt.plot(x_list, p(x_list), "r--")

    plt.title(prop2_str
              + " vs "
              + prop1_str
              + " for "
              + atomName
              + " scatter plot")
    plt.xlabel(prop1_str)
    plt.ylabel(prop2_str)

    plt.show()

def plot_prop_different_elements(prop1_str,
                                 prop2_str,
                                 struct_str,
                                 database,
                                 integrator,
                                 potential):
    file = open("sim_properties.csv", "r")
    x_list = []
    y_list = []
    reader = csv.DictReader(file, delimiter = ";")
    sort = sorted(reader, key = operator.itemgetter("Element"))
    element = sort[0]["Element"]
    row = 0
    column = 0
    simulation_matrix = np.empty([1000,1000])
    element_list = ["Symbols"]
    for listrow in sort:
        if listrow["Database"] == database and listrow["Potential"] == potential:
            if listrow["Integrator"] == integrator and listrow["Struct"].startswith(struct_str):
                if listrow["Element"] == element:
                    if column == 0:
                        element_list.append(listrow["Element"])
                        element_list.append("Empty")
                        simulation_matrix[row,column] = listrow[prop1_str]
                        simulation_matrix[row + 1, column] = listrow[prop2_str]
                        column += 1
                    else:
                        simulation_matrix[row,column] = listrow[prop1_str]
                        simulation_matrix[row + 1, column] = listrow[prop2_str]
                        column +=1
                        element = listrow["Element"]
                else:
                    row += 2
                    column = 0
                    element = listrow["Element"]

    i = 0
    n = 0
    while i < len(simulation_matrix):
        if simulation_matrix[i,0] != 0:
            x_list = list(filter(lambda a: a != 0, simulation_matrix[i]))
            y_list = list(filter(lambda a: a != 0, simulation_matrix[i+1]))
            plt.scatter(x_list, y_list, label = element_list[i+1])
            #Add trendline to scatterplots.
            #z = np.polyfit(x_list, y_list, 1)
            #p = np.poly1d(z)
            #plt.plot(x_list, p(x_list), "r--")
            i += 2
        else:
            break

    plt.title(prop2_str
              + " vs "
              + prop1_str
              + " for all simulated elements"
              + " scatter plot")
    plt.xlabel(prop1_str)
    plt.ylabel(prop2_str)
    plt.legend()
    plt.show()
    

def hist_prop_per_simulation(atomName,
                             prop_str,
                             struct_str,
                             database,
                             integrator,
                             potential):
    file = open("sim_properties.csv", "r")
    reader = csv.DictReader(file, delimiter=";")
    count = 1
    plot_list = []
    for row in reader:
        if row["Database"] == database and row["Potential"] == potential:
            if row["Integrator"] == integrator and row["Struct"].startswith(struct_str):
                if row["Element"] == atomName:
                    plot_list.append(float(row[prop_str]))
                    count += 1
        else:
            return(None)
    plt.hist(plot_list, bins = 10)

    plt.title(prop_str + " histogram")
    plt.xlabel(prop_str)
    
    plt.show()

"""
Available databases: "Materials Project" or "ASE"
Availablable integrators: "Velocity-Verlet" or "Langevin"
Available potentials: "EMT", "Leonnard-Jones" or  specific key for chosen kim potential as a string.
Available elements: Check sim_properties.csv columnt "Elements" to see which elements have been simulated.
For EMT: Cu, Ni, Au, Ag, Al, Pd, Pt (H, C, N, O only for fun appearently).
Available properties for plotting: "MSD", "D", "L", "SHC", "Int_T", "E_coh", "Int_P", "Bulk_mod".
Structures: Look in the sim_properties.csv too see which structures are saved.
"""

#Run plot functions. Will appear one after another.

plot_prop_vs_time("properties.csv", "MSD")
plot_prop_per_simulation("Cu", "CUB", "Int_T", "L", "ASE", "Velocity-Verlet", "EMT")
plot_prop_different_elements("Int_T", "E_coh", "CUB", "ASE", "Velocity-Verlet", "EMT")
hist_prop_per_simulation("Cu", "MSD", "CUB", "ASE", "Velocity-Verlet", "EMT")
