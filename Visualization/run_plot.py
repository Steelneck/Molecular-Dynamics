import csv
import matplotlib.pyplot as plt

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
   
    plt.title(prop_str + " vs time scatter plot")
    plt.xlabel("Time [fs]")
    plt.ylabel(prop_str)

    plt.show()

def plot_prop_per_simulation(AtomName, prop1_str, prop2_str):
    file = open("properties_" + AtomName + ".csv", "r")
    reader = csv.DictReader(file, delimiter=";")
    count = 1
    x_list = []
    y_list = []
    for row in reader:
        x_list.append(float(row[prop1_str]))
        y_list.append(float(row[prop2_str]))
        count += 1

    plt.scatter(x_list, y_list)

    plt.title(prop2_str
              + " vs "
              + prop1_str
              + " for "
              + AtomName
              + " scatter plot")
    plt.xlabel(prop1_str)
    plt.ylabel(prop2_str)

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

plot_prop_vs_time("properties.csv", "MSD")
#plot_prop_per_simulation("Cu", "Int_P", "Int_T")
#hist_prop_per_simulation("Cu", "MSD")


