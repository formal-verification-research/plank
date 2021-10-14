# Description
'''
The purpose of Anastomosis is to kill the EC if they run into a different cell line
'''


# Function
def anastomosis(anastomotic, yPosition, workspace, fileDeaths, deathTime, occupied, cell, currentTimeStep, xPosition, cellLine):

    if anastomotic is True:
        if yPosition[cell][currentTimeStep + 1] > 0:
            if workspace[yPosition[cell][currentTimeStep + 1]][xPosition[cell][currentTimeStep + 1]] != 0 and \
                    workspace[yPosition[cell][currentTimeStep + 1]][xPosition[cell][currentTimeStep + 1]] != cellLine[cell]:

                fileDeaths.write("\n\nTime: " + str(currentTimeStep) + "\n")
                fileDeaths.write("Cell: " + str(cell) + "\n")
                fileDeaths.write("Cell ran into another capillary" + "\n")
                deathTime[cell] = currentTimeStep + 1
                occupied[yPosition[cell][currentTimeStep + 1]][xPosition[cell][currentTimeStep + 1]] -= 1

    workspace[yPosition[cell][currentTimeStep]][xPosition[cell][currentTimeStep]] = cellLine[cell]
    workspace[yPosition[cell][currentTimeStep + 1]][xPosition[cell][currentTimeStep + 1]] = 100

    return
