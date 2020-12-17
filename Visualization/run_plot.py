import csv, operator, sys, os
import matplotlib.pyplot as plt
import numpy as np

def plot_prop_vs_time(csvFileName, prop_str):

    """
    Function which plots the time evolution of property specified by input prop_str. csvFileName is the .csv-file where the time-evolution of the latest run simulation is saved.
    """
    try:
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
        
        #plt.scatter(x_list,y_list)

        #Adds trendline to the plot.
       # z = np.polyfit(x_list, y_list, 1)
       # p = np.poly1d(z)
       # plt.plot(x_list, p(x_list), "r--")

       # plt.title(prop_str + " vs time scatter plot")
       # plt.xlabel("Time [fs]")
       # plt.ylabel(prop_str)
        
        #Set interval for the y-axis to not have the plot zoomed in too much.
       # plt.ylim(0,0.1)
        
        #plt.show()
        return(x_list)
        
    except Exception as e:
        print("An error occured when plotting a property vs time:")
        exc_type, exc_obj, exc_traceBack = sys.exc_info()
        fname = os.path.split(exc_traceBack.tb_frame.f_code.co_filename)[1]
        print("Error type:", exc_type, "; Message:", e, "; In file:", fname, "; On line:", exc_traceBack.tb_lineno)
        return(None)

def plot_prop_per_simulation(atomName,
                             struct_str,
                             prop1_str,
                             prop2_str,
                             database,
                             integrator,
                             potential):
    """
    Functions that creates a scatterplot with prop_1 str on the x-axis and prop2_str on the y-axis for element with name atomName taken from database, simulated with integrator and potential. Not sure if this one will be necessary.
    """

    try:
        file = open("sim_properties.csv", "r")
        reader = csv.DictReader(file, delimiter=";") #file is of dictionary format so the reader needs to be a dictreader. ";" separates the columns.
        count = 1
        x_list = []
        y_list = []
        #loop over all elements in the dictionary.
        for row in reader:
        #Check if the current row has the correct simulation input.
            if row["Database"] == database and row["Integrator"] == integrator:
                if row["Potential"] == potential:
                    if row["Element"] == atomName and row["Struct"].startswith(struct_str):
                    #If all criteria is fullfilled, the value on the current row in column prop1_str and prop2_str is appended to two separate lists.
                        x_list.append(float(row[prop1_str]))
                        y_list.append(float(row[prop2_str]))
                        count += 1             
        #plt.scatter(x_list, y_list)

        #These lines of code are only needed if a trendline is requested.
        #z = np.polyfit(x_list, y_list, 1)
        #p = np.poly1d(z)
        #plt.plot(x_list, p(x_list), "r--")

        #plt.title(prop2_str
        #          + " vs "
        #          + prop1_str
        #          + " for "
        #          + atomName
        #          + " scatter plot")
        #plt.xlabel(prop1_str)
        #plt.ylabel(prop2_str)
        
        #plt.show()
        return(x_list)

    except Exception as e:
        print("An error occured when plotting properties for all simulations:")
        exc_type, exc_obj, exc_traceBack = sys.exc_info()
        fname = os.path.split(exc_traceBack.tb_frame.f_code.co_filename)[1]
        print("Error type:", exc_type, "; Message:", e, "; In file:", fname, "; On line:", exc_traceBack.tb_lineno)
        return(None)
    
def plot_prop_different_elements(prop1_str, prop2_str, struct_str,
                                 database, integrator):
    """
    Function which plots two properties (specified by prop1_str and prop2_str taken from .csv-file sim_properties.csv) and creates a scatter plot of all atom-types and potentials separating them with colors and markers. Currently very ugly function but lack of time stops me from cleaning it up. It does the job.
    """

    try:
        file = open("sim_properties.csv", "r")
        marker_list = ["+", "o", "*", "x", "v", "d", "^", "s", "2"] #Make a list with markers that can be looped through when plotting multiple datatypes in the same scatterplot.
        x_list = []
        y_list = []
        reader = csv.DictReader(file, delimiter = ";")
        sort = sorted(reader, key = operator.itemgetter("Element", "Potential")) #Create sorted list after element and then potential to make it easier to loop through
        element = sort[0]["Element"] #Save first element to compare.
        potential = sort[0]["Potential"] #Save first potential to compare.
        row = 0 #Row index for simulation_matrix
        column = 0 #Column index for simulation_matrix
        simulation_matrix = np.empty([1000,5000]) #Make a list that is big enough to handle the largest datasets that we will run. This will be trimmed depending on how many values we append.
        element_list = ["Symbols followed by potential"] #List that will contain the symbols of each element followed by the potential both as strings.
        for listrow in sort:
            #Check criteria for correct input.
            if listrow["Database"] == database and listrow["Integrator"]== integrator:
                if listrow["Struct"].startswith(struct_str):
                    #This is where the current config for "Potential and "Element" are compared to the previous iteration to see if we need too switch rows.
                    if listrow["Element"] == element:
                        if listrow["Potential"] == potential:
                            #If column is 0 we want to save dict values under "Element" and "Potential" to element_list and also append the values to simulation matrix.
                            if column == 0:
                                element_list.append(listrow["Element"])
                                element_list.append(listrow["Potential"])
                                simulation_matrix[row,column] = listrow[prop1_str]
                                simulation_matrix[row + 1, column] = listrow[prop2_str]
                                column +=1
                            else:
                                #If column != 0 we only want to add values to simulation_matrix.
                                simulation_matrix[row,column] = listrow[prop1_str]
                                simulation_matrix[row + 1, column] = listrow[prop2_str]
                                column +=1
                        #If the current row has a new potential, we want to jump down in simulation_matrix and then save the values to that row on column index 0.
                        else:
                            row += 2
                            column = 0
                            element_list.append(listrow["Element"])
                            element_list.append(listrow["Potential"])
                            simulation_matrix[row,column] = listrow[prop1_str]
                            simulation_matrix[row + 1, column] = listrow[prop2_str]
                            column +=1
                    #If the current row has a new element we also want to jump down in simulation_matrix.
                    else:
                        row += 2
                        column = 0
                        element_list.append(listrow["Element"])
                        element_list.append(listrow["Potential"])
                        simulation_matrix[row,column] = listrow[prop1_str]
                        simulation_matrix[row + 1, column] = listrow[prop2_str]
                        column+=1
            #Save the current values to compare in next iteration.
            potential = listrow["Potential"]
            element = listrow["Element"]
        
        i = 0
        marker_counter = 0 #Used to change markers in marker_list for different plots.
        #Now loop through all rows of simulation_matrix and make separate plots for each element and potential, separated by color and marker type.
        while i < len(simulation_matrix):
            if simulation_matrix[i,0] != 0 and i+2<len(element_list):
                x_list = list(filter(lambda a: a != 0, simulation_matrix[i]))
                y_list = list(filter(lambda a: a != 0, simulation_matrix[i+1]))
                #Make sure the label says what element and potential the marker and color represents.
                plt.scatter(x_list, y_list, label = element_list[i+1] + " with potential " + element_list[i+2], marker = marker_list[marker_counter])
                #Add trendline to scatterplots.
                #z = np.polyfit(x_list, y_list, 1)
                #p = np.poly1d(z)
                #plt.plot(x_list, p(x_list), "r--")
                i += 2
                marker_counter += 1
            else:
                break

        plt.title(prop2_str
                 + " vs "
                  + prop1_str
                  + " for all simulated elements"
                  + " scatter plot")
        plt.xlabel(prop1_str)
        plt.ylabel(prop2_str)
        plt.ylim(0,0.1)
        plt.legend()
        plt.show()
        print(x_list)
        return(x_list)
    
    
    except Exception as e:
        print("An error occured when plotting properties for all elements and potentials:")
        exc_type, exc_obj, exc_traceBack = sys.exc_info()
        fname = os.path.split(exc_traceBack.tb_frame.f_code.co_filename)[1]
        print("Error type:", exc_type, "; Message:", e, "; In file:", fname, "; On line:", exc_traceBack.tb_lineno)
        return(None)

    
def hist_prop_per_simulation(atomName,
                             prop_str,
                             struct_str,
                             database,
                             integrator,
                             potential):
    """
    Function that makes a histogram plot of a chosen property (prop_str) for atom with atomName and structure struct_str, simulated with database, integrator and potential.
    """
    try:
        file = open("sim_properties.csv", "r")
        reader = csv.DictReader(file, delimiter=";")
        count = 1
        plot_list = []
        #Similar to previous function. Loop through sim_properties with dictionary reader and add property to plot_list if input criteria is fulfilled.
        for row in reader:
            if row["Database"] == database and row["Potential"] == potential:
                if row["Integrator"] == integrator and row["Struct"].startswith(struct_str):
                    if row["Element"] == atomName:
                        plot_list.append(float(row[prop_str]))
                        count += 1
                        
        plt.hist(plot_list, bins = 10)

        plt.title(prop_str + " histogram")
        plt.xlabel(prop_str)
        plt.show()
        return(plot_list)
        
        
    except Exception as e:
        print("An error occured when plotting histogram of property:")
        exc_type, exc_obj, exc_traceBack = sys.exc_info()
        fname = os.path.split(exc_traceBack.tb_frame.f_code.co_filename)[1]
        print("Error type:", exc_type, "; Message:", e, "; In file:", fname, "; On line:", exc_traceBack.tb_lineno)
        return(None)

"""
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

#plot_prop_vs_time("properties.csv", "MSD")
#plot_prop_per_simulation("Cu", "CUB", "Int_T", "L", "ASE", "Velocity-Verlet", "EMT")
#plot_prop_different_elements("Int_T", "MSD", "CUB", "ASE", "Velocity-Verlet")
#hist_prop_per_simulation("Cu", "MSD", "CUB", "ASE", "Velocity-Verlet", "EMT")
