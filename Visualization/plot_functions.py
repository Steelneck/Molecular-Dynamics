import csv, operator, sys, os, json
import matplotlib.pyplot as plt
import numpy as np

def plot_prop_vs_time(jsonFileName):
    
    """
    Function which plots the time evolution of property specified by input prop_str. csvFileName is the .csv-file where the time-evolution of the latest run simulation is saved.
    """
    
    try:
        #Load relevant data from input file "plot_input.json"
        with open(jsonFileName, "r") as json_file:
            data = json.load(json_file)
            for row in data["Scatter plot, time evolution of single property"]:
                prop_str = row["Property"]
                plot_check = row["Plot check"]
                #The following makes it possible to scale the plot x- & y-axis
                scale_axis_check = row["Scale axis check"]
                y_min = row["y-min"]
                y_max = row["y-max"]
                x_min = row["x-min"]
                x_max = row["x-max"]
                
        file = open("properties.csv", "r")
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
        
        plt.title(prop_str + " vs time scatter plot")
        plt.xlabel("Time [fs]")
        plt.ylabel(prop_str)
        
        #Set interval for the y-axis to not have the plot zoomed in too much.
        if scale_axis_check:
            if x_min != x_max:
                plt.xlim(x_min,x_max)
            if y_min != y_max:
                plt.ylim(y_min,y_max)
        
        if plot_check:
            plt.show()
        plt.close()
        return(x_list)
        
    except Exception as e:
        print("An error occured when plotting a property vs time:")
        exc_type, exc_obj, exc_traceBack = sys.exc_info()
        fname = os.path.split(exc_traceBack.tb_frame.f_code.co_filename)[1]
        print("Error type:", exc_type, "; Message:", e, "; In file:", fname, "; On line:", exc_traceBack.tb_lineno)
        return(None)

def plot_prop_per_simulation(jsonFileName):
    """
    Functions that creates a scatterplot with prop_1 str on the x-axis and prop2_str on the y-axis for element with name atomName taken from database, simulated with integrator and potential. Not sure if this one will be necessary.

    atomName, struct_str, prop1_str, prop2_str, database, integrator, potential
    """

    try:
        #Load relevant data from input file "plot_input.json"
        with open(jsonFileName, "r") as json_file:
            data = json.load(json_file)
            for row in data["Scatter plot, two properties for a specific element"]:
                prop1_str = row["Property x-axis"]
                prop2_str = row["Property y-axis"]
                atomName = row["Element"]
                struct_str = row["Structure"]
                database = row["Database"]
                integrator = row["Integrator"]
                potential = row["Potential"]
                plot_check = row["Plot check"]
                #The following makes it possible to scale the plot x- & y-axis
                scale_axis_check = row["Scale axis check"]
                y_min = row["y-min"]
                y_max = row["y-max"]
                x_min = row["x-min"]
                x_max = row["x-max"]
                
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
        if scale_axis_check:
            if x_min != x_max:
                plt.xlim(x_min,x_max)
            if y_min != y_max:
                plt.ylim(y_min,y_max)
        if plot_check:
            plt.show()
        plt.close()
        return(x_list)

    except Exception as e:
        print("An error occured when plotting properties for all simulations:")
        exc_type, exc_obj, exc_traceBack = sys.exc_info()
        fname = os.path.split(exc_traceBack.tb_frame.f_code.co_filename)[1]
        print("Error type:", exc_type, "; Message:", e, "; In file:", fname, "; On line:", exc_traceBack.tb_lineno)
        return(None)
    
def plot_prop_all_elements(jsonFileName):
    """
    Function which plots two properties (specified by prop1_str and prop2_str taken from .csv-file sim_properties.csv) and creates a scatter plot of all atom-types and potentials separating them with colors and markers. Currently very ugly function but lack of time stops me from cleaning it up. It does the job.
    prop1_str, prop2_str, struct_str, database, integrator
    """

    try:
        #Load relevant data from input file "plot_input.json"
        with open(jsonFileName, "r") as json_file:
            data = json.load(json_file)
            for row in data["Scatter plot, two properties for all elements"]:
                prop1_str = row["Property x-axis"]
                prop2_str = row["Property y-axis"]
                struct_str = row["Structure"]
                database = row["Database"]
                integrator = row["Integrator"]
                plot_check = row["Plot check"]
                #The following makes it possible to scale the plot x- & y-axis
                scale_axis_check = row["Scale axis check"]
                y_min = row["y-min"]
                y_max = row["y-max"]
                x_min = row["x-min"]
                x_max = row["x-max"]
            
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
        if scale_axis_check:
            if x_min != x_max:
                plt.xlim(x_min,x_max)
            if y_min != y_max:
                plt.ylim(y_min,y_max)
        plt.legend()
        if plot_check:
            plt.show()
        plt.close()
        return(x_list)
    
    
    except Exception as e:
        print("An error occured when plotting properties for all elements and potentials:")
        exc_type, exc_obj, exc_traceBack = sys.exc_info()
        fname = os.path.split(exc_traceBack.tb_frame.f_code.co_filename)[1]
        print("Error type:", exc_type, "; Message:", e, "; In file:", fname, "; On line:", exc_traceBack.tb_lineno)
        return(None)

    
def hist_prop_per_simulation(jsonFileName):
    """
    Function that makes a histogram plot of a chosen property (prop_str) for atom with atomName and structure struct_str, simulated with database, integrator and potential.
    atomName, prop_str, struct_str, database, integrator, potential
    """
    try:
        #Load relevant data from input file "plot_input.json"
        with open(jsonFileName, "r") as json_file:
            data = json.load(json_file)
            for row in data["Histogram, single property from all simulations"]:
                atomName = row["Element"]
                prop_str = row["Property"]
                struct_str = row["Structure"]
                database = row["Database"]
                integrator = row["Integrator"]
                potential = row["Potential"]
                plot_check = row["Plot check"]
            
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
        if plot_check:
            plt.show()
        plt.close()
        return(plot_list)
        
        
    except Exception as e:
        print("An error occured when plotting histogram of property:")
        exc_type, exc_obj, exc_traceBack = sys.exc_info()
        fname = os.path.split(exc_traceBack.tb_frame.f_code.co_filename)[1]
        print("Error type:", exc_type, "; Message:", e, "; In file:", fname, "; On line:", exc_traceBack.tb_lineno)
        return(None)
