import csv, operator
import matplotlib.pyplot as plt
import numpy as np

def plot_prop_vs_time(csvFileName, prop_str):
    file = open(csvFileName, "r")
    reader = csv.DictReader(file, delimiter=";")
    count = 1
    x_list = []
    y_list = []
    for row in reader:
        x_list.append(float(row["Time"]))
        y_list.append(float(row[prop_str]))
        count += 1
        
    plt.scatter(x_list,y_list)

    z = np.polyfit(x_list, y_list, 1)
    p = np.poly1d(z)
    plt.plot(x_list, p(x_list), "r--")

    plt.title(prop_str + " vs time scatter plot")
    plt.xlabel("Time [fs]")
    plt.ylabel(prop_str)
    plt.ylim(0,0.1)

    plt.show()

def plot_prop_per_simulation(atomName,
                             struct,
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
            if ["Potential"] == potential:
                if row["Element"] == atomName and row["Struct"].startswith(struct):
                    x_list.append(float(row[prop1_str]))
                    y_list.append(float(row[prop2_str]))
                    count += 1
                    
    plt.scatter(x_list, y_list)

    z = np.polyfit(x_list, y_list, 1)
    p = np.poly1d(z)
    plt.plot(x_list, p(x_list), "r--")

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
            #z = np.polyfit(x_list, y_list, 1)
            #p = np.poly1d(z)
            #plt.plot(x_list, p(x_list), "r--")
            i += 2
        else:
            break
    plt.legend()
    plt.show()
    

def hist_prop_per_simulation(AtomName, prop_str):
    file = open("properties_" + AtomName + ".csv", "r")
    reader = csv.DictReader(file, delimiter=";")
    count = 1
    plot_list = []
    for row in reader:
        plot_list.append(float(row[prop_str]))
        count += 1
        plt.hist(plot_list, bins = 10)

    plt.title(prop_str + " histogram")
    plt.xlabel(prop_str)
    
    plt.show()

"""
Available databases: "Materials Project" or "ASE"
Availablable integrators: "Velocity-Verlet" or "Langevin"
Available potentials: "EMT, "Leonnard-Jones" or  specific key for chosen kim potential as a string.
Available elements: "Cu", 
"""

plot_prop_vs_time("properties.csv", "MSD")
#plot_prop_per_simulation("Cu", "FCC", "Int_T", "L")
plot_prop_different_elements("Int_T", "E_coh", "CUB", "ASE", "Velocity-Verlet", "EMT")
#hist_prop_per_simulation("Cu", "MSD")
