import csv
import matplotlib.pyplot as plt

def plot_from_csvFile(csvFileName, x_str, y_str):
    file = open(csvFileName, "r")
    reader = csv.DictReader(file, delimiter=";")
    count = 1
    x_list = []
    y_list = []
    for row in reader:
        x_list.append(float(row[x_str]))
        y_list.append(float(row[y_str]))
        count += 1
        
    plt.scatter(x_list,y_list)
   
    plt.title("Self-diffusion vs MSD scatter plot")
    plt.xlabel(x_str)
    plt.ylabel(y_str)

    plt.show()

plot_from_csvFile("properties.csv", "Time", "MSD")
