# Description
'''
The purpose of Division is to divide the cells in a cross pattern like is described in the plank paper
'''


# Function
def division(x, y, xPosition, yPosition, occupied, deathTime, birthTime, divideTime, numberOfCells, cell, currentTimeStep,
             totalNumberOfTimeSteps, divideProbability, fileDivisions, cellLine):

    # Divide the cell and create a new one
    xPosition[numberOfCells][currentTimeStep + 1] = x
    yPosition[numberOfCells][currentTimeStep + 1] = y
    occupied[y][x] += 1
    cellLine[numberOfCells] = cellLine[cell]
    deathTime[numberOfCells] = totalNumberOfTimeSteps - 1
    birthTime[numberOfCells] = currentTimeStep
    divideTime[cell] = currentTimeStep
    divideTime[numberOfCells] = currentTimeStep
    numberOfCells += 1

    # Get some data into the divison file
    fileDivisions.write("\n\nTime: " + str(currentTimeStep) + "\n")
    fileDivisions.write("Cell: " + str(cell) + "\n")
    fileDivisions.write("Probability of Division: " + str(divideProbability) + "\n")

    return
