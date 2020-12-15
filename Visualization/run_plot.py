import csv
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

    plt.show()

def plot_prop_per_simulation(atomName, struct, prop1_str, prop2_str):
    file = open("sim_properties.csv", "r")
    reader = csv.DictReader(file, delimiter=";")
    count = 1
    x_list = []
    y_list = []
    for row in reader:
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

def plot_prop_different_elements(struct_str, temp_str,  prop1_str, prop2_str):
    file = open("sim_properties.csv", "r")
    reader = csv.DictReader(file, delimiter = ";")
    count = 1
    x_list = []
    y_list = []
    element_list = []
    for row in reader:
        if row["Init_T"] == temp_str and row["Struct"].startswith(struct_str):
            x_list.append(float(row[prop1_str]))
            y_list.append(float(row[prop2_str]))
            element_list.append(row["Element"])
            count += 1
        else:
            count += 1
    
    plt.scatter(x_list, y_list)
    for i, txt in enumerate(element_list):
        plt.annotate(txt, (x_list[i], y_list[i]))
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

#plot_prop_vs_time("properties.csv", "MSD")
#plot_prop_per_simulation("Cu", "FCC", "Int_T", "L")
plot_prop_different_elements("FCC", "300", "Int_T", "L")
#hist_prop_per_simulation("Cu", "MSD")
