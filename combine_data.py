# this file combines all the collected data from each run and makes excel sheets and text files

import pandas as pd
import glob
import sys
import numpy as np
from parameter_vault import NUMBER_OF_RUNS


def combineData():
    # this function is not used, but contained important information of path file
    file = r"C:\Users\cassa\Desktop\Results\test"
    files = (glob.glob(file+"/**/Density"))
    number = 0
    for file in files:
        file = pd.read_csv(file, delim_whitespace=True)
        f = pd.DataFrame(file)
        with pd.ExcelWriter(r"C:\Users\cassa\Desktop\Results\test\test.xlsx", mode="a", engine="openpyxl") as writer:
            f.to_excel(writer, sheet_name=str(number))
        number += 1


def combineDensity():
    # open each file and pull the needed information
    file = r"C:\Users\cassa\Desktop\Results\test"
    files = (glob.glob(file + "/**/Density"))
    allData = {}
    fileIndex = 0
    for fl in files:
        openFile = open(fl, "r")
        for line in openFile.readlines():
            line = line.split("\t")
            if len(line) != 2:
                print(f"Line {line} is misformatted! Ignoring line.", file=sys.stderr)
                continue
            time = line[0]
            density = line[1]
            if not time in allData:
                allData[time] = ["" for i in range(fileIndex + 1)]
                allData[time][fileIndex] = density
            else:
                currentRow = allData[time]
                while len(currentRow) <= fileIndex:
                    currentRow.append("")
                currentRow[fileIndex] = density
        fileIndex += 1
        openFile.close()

    # create headers and titles of each column based on the file path
    headers = ["Time (Hours)"]
    for f in files:
        headers.append(str(f))
    headers.append('average_density')
    data = [headers]

    # initialize some vectors for gathering information
    endTimes = []
    endDensity = []
    skipValue = []
    skipVal = False
    ave_density = []

    # find the ending density in each simulation
    for time, densityRow in allData.items():
        currentRow = [time]
        for value in currentRow:
            value = float(value)
            if not value % 0.5 == 0:
                endTimes.append(value)
                skipVal = True
                for density in densityRow:
                    if not density == "":
                        endDensity.append(density)
            skipValue.append(skipVal)
            skipVal = False

    # create another column and add in the average density
    for time, densityRow in allData.items():
        currentRow=[time]
        for density in densityRow:
            for i in range(len(currentRow)):
                if density == '':
                    if len(endDensity) == 0:
                        pass
                    else:
                        density = endDensity[len(currentRow)-1]
            currentRow.append(density)
            if len(currentRow) == NUMBER_OF_RUNS + 1:
                density_num = 0
                for i in currentRow:
                    if i is currentRow[0]:
                        pass
                    else:
                        if i == '':
                            pass
                        else:
                            i = float(i)
                            density_num = (density_num + i)
                currentRow.append(float(density_num/NUMBER_OF_RUNS))
                ave_density.append(float(density_num/NUMBER_OF_RUNS))

    # check for the length of each row and fill in any blank spaces
        while len(currentRow) < NUMBER_OF_RUNS + 2:
            if len(currentRow) == NUMBER_OF_RUNS + 1:
                currentRow.append(ave_density[-1])
            else:
                currentRow.append(endDensity[len(currentRow)-1])
        if float(currentRow[0]) % 0.5 == 0 or float(currentRow[0]) == max(endTimes):
            data.append(currentRow)

    # recalculate the average density
    for i in range(1, len(data)):
        sum = 0
        for j in range(1, NUMBER_OF_RUNS+1):
            data[i][j] = float(data[i][j])
            sum += data[i][j]
        data[i][-1] = sum / NUMBER_OF_RUNS

    # calculate the change in density and append it to the end of the data table
    change_density_array = []
    change_density_array.append("AverageChangeinDensity")
    change_density_array.append(0)
    for j in range(2, len(data)):
        if data[j][-1] == str:
            change_density = "average change in density"
        elif data[j-1][-1] == str:
            change_density = "average change in density"
        else:
            change_density = float(data[j][NUMBER_OF_RUNS + 1]) - float(data[j-1][NUMBER_OF_RUNS + 1])
            change_density_array.append(float(change_density))
    for i in range(len(change_density_array)):
        data[i].append(change_density_array[i])

    endDensity2 = []
    for i in endDensity:
        i = float(i)
        endDensity2.append(i)
    average_end_density = np.mean(endDensity2)
    average_end_time = np.mean(endTimes)
    max_end_time = np.max(endTimes)
    file_average = open("Average Density", 'w')
    file_average.write("Average End Density (%):  " + str(average_end_density) + "\n" + "Average End Time (hrs):  " +
                                str(average_end_time) + "\n" + "Max End Time(hrs): " + str(max_end_time))

    frame = pd.DataFrame(data, index=None, columns=headers)

    print(frame.tail(5))
    with pd.ExcelWriter(r"C:\Users\cassa\Desktop\Results\test\test_other.xlsx", mode="w", engine="openpyxl") as writer:
        frame.to_excel(writer, sheet_name="Aggregate Data")

    # collect the average density and the average change over time, to use in the combineAll function
    ave_density = []
    change_density_array = []
    for i in range(1, len(data)):
        ave_density.append(data[i][-2])
        change_density_array.append(data[i][-1])

    return ave_density, change_density_array


def combineNodes():
    file2 = r"C:\Users\cassa\Desktop\Results\test"
    files2 = (glob.glob(file2 + "/**/Nodes"))
    allData2 = {}
    fileIndex2 = 0

    # open all files and read through nodes information
    for fl in files2:
        openFile2 = open(fl, "r")
        for line in openFile2.readlines():
            line = line.split("\t")
            if len(line) != 2:
                print(f"Line {line} is misformatted! Ignoring line.", file=sys.stderr)
                continue
            time = line[0]
            nodes = line[1]
            if not time in allData2:
                allData2[time] = ["" for i in range(fileIndex2 + 1)]
                allData2[time][fileIndex2] = nodes
            else:
                currentRow = allData2[time]
                while len(currentRow) <= fileIndex2:
                    currentRow.append("")
                currentRow[fileIndex2] = nodes
        fileIndex2 += 1
        openFile2.close()
    headers2 = ["Time (Hours)"]
    for f in files2:
        headers2.append(str(f))
    headers2.append('AverageNodes')
    data_nodes = [headers2]

    # initiate vectors for gathering information
    endTimesNodes = []
    endNodes = []
    skipValueNodes = []
    skipValNodes = False
    ave_nodes = []

    # find the final amount of nodes at the end of each simulation
    for time, nodesRow in allData2.items():
        currentRowNodes = [time]
        for value in currentRowNodes:
            value = float(value)
            if not value % 0.5 == 0:
                endTimesNodes.append(value)
                skipValNodes = True
                for nodes in nodesRow:
                    if not value % 0.5 == 0 and not nodes == "":
                        endNodes.append(nodes)
            skipValueNodes.append(skipValNodes)
            skipValNodes = False

    # fill in any blank spaces from simulations that ended early and add in another row for average nodes
    for time, nodesRow in allData2.items():
        currentRowNodes = [time]
        for density in nodesRow:
            for i in range(len(currentRowNodes)):
                if density == '':
                    density = endNodes[len(currentRowNodes)-1]
            currentRowNodes.append(density)
            if len(currentRowNodes) == NUMBER_OF_RUNS + 1:
                nodes_num = 0
                for i in currentRowNodes:
                    if i is currentRowNodes[0]:
                        pass
                    else:
                        i = float(i)
                        nodes_num = (nodes_num + i)
                currentRowNodes.append(nodes_num/NUMBER_OF_RUNS)
                ave_nodes.append(nodes_num/NUMBER_OF_RUNS)

    # only add the rows to the end table if they have the right length, fill in any blanks with the appropriate number
        while len(currentRowNodes) < NUMBER_OF_RUNS + 2:
            if len(currentRowNodes) == NUMBER_OF_RUNS + 1:
                currentRowNodes.append(ave_nodes[-1])
            else:
                currentRowNodes.append(endNodes[len(currentRowNodes)-1])
        if float(currentRowNodes[0]) % 0.5 == 0 or float(currentRowNodes[0]) == max(endTimesNodes):
            data_nodes.append(currentRowNodes)

    # recalculate the correct average nodes
    for i in range(1, len(data_nodes)):
        sum = 0
        for j in range(1, NUMBER_OF_RUNS+1):
            data_nodes[i][j] = int(data_nodes[i][j])
            sum += data_nodes[i][j]
        data_nodes[i][-1] = sum / NUMBER_OF_RUNS
        ave_nodes.append(data_nodes[i][-1])

    # calculate the change in nodes and append it to the end of data
    change_nodes_array = []
    change_nodes_array.append("AverageChangeinNodes")
    change_nodes_array.append(0)
    for j in range(2, len(data_nodes)):
        if data_nodes[j][-1] == str:
            change_nodes = "average change in nodes"
        elif data_nodes[j-1][-1] == str:
            change_nodes = "average change in nodes"
        else:
            change_nodes = float(data_nodes[j][NUMBER_OF_RUNS + 1]) - float(data_nodes[j-1][NUMBER_OF_RUNS + 1])
            change_nodes_array.append(float(change_nodes))
    for i in range(len(change_nodes_array)):
        data_nodes[i].append(change_nodes_array[i])

    # collect average end times and max run time
    endNodes2 = []
    for i in endNodes:
        i = float(i)
        endNodes2.append(i)
    average_end_nodes = np.mean(endNodes2)
    average_end_time = np.mean(endTimesNodes)
    max_end_time_node = max(endTimesNodes)
    file_average = open("Average Nodes", 'w')
    file_average.write("Average End Nodes (%):  " + str(average_end_nodes) + "\n" + "Average End Time (hrs):  " +
                                str(average_end_time) + "\n" + "Max End Time: " + str(max_end_time_node))

    # convert data table to excel format for easy reading with pandas
    frame = pd.DataFrame(data_nodes, index=None, columns=headers2)
    print(frame)
    with pd.ExcelWriter(r"C:\Users\cassa\Desktop\Results\test\test_nodes.xlsx", mode="w", engine="openpyxl") as writer:
        frame.to_excel(writer, sheet_name="Aggregate Data")

    # collect the average number of nodes and the average change over time, to use in the combineAll function
    ave_nodes = []
    change_nodes_array = []
    for i in range(1, len(data_nodes)):
        ave_nodes.append(data_nodes[i][-2])
        change_nodes_array.append(data_nodes[i][-1])

    return ave_nodes, change_nodes_array


def combineAll():

    # run the combineDensity and combineNodes functions and combine averages from both onto a separate sheet
    ave_density, change_density_array = combineDensity()
    ave_nodes, change_nodes_array = combineNodes()

    density_nodes_array = [([float(0)]*4) * len(ave_density)]
    density_nodes_array = np.reshape(density_nodes_array, (len(ave_density), 4))

    for i in range(len(ave_density)):
        density_nodes_array[i][0] = ave_density[i]
    for i in range(len(ave_density)):
        density_nodes_array[i][1] = change_density_array[i]
    for i in range(len(ave_density)):
        density_nodes_array[i][2] = ave_nodes[i]
    for i in range(len(ave_density)):
        density_nodes_array[i][3] = change_nodes_array[i]

    headers3= ["Average Density over Time", "Change in Density", "Average Nodes over Time", "Change in Nodes"]
    density_nodes = pd.DataFrame(density_nodes_array, index=None, columns=headers3)

    print(density_nodes)

    with pd.ExcelWriter(r"C:\Users\cassa\Desktop\Results\test\test_nodesanddensity.xlsx", mode="w", engine="openpyxl") \
            as writer:
        density_nodes.to_excel(writer, sheet_name="Aggregate Data")

if __name__=="__main__":
    combineAll()
    # combineDensity()

