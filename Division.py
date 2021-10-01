# Description
'''
The purpose of Division is to divide the cells in a cross pattern like is described in the plank paper
'''


# Function
def division(x, y, xPosition, yPosition, occupied, deathTime, birthTime, divideTime, numberOfCells, cell, currentTimeStep,
             totalNumberOfTimeSteps, divideProbability, file2, cellLine):

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
    file2.write("\n\nTime: " + str(currentTimeStep) + "\n")
    file2.write("Cell: " + str(cell) + "\n")
    file2.write("Probability of Division: " + str(divideProbability) + "\n")

    return
