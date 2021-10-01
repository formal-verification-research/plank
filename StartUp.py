# Description
'''
The StartUp file is to create the fields and to place the cells inside the simulation
'''


# Function
def startUp(xSteps, numberOfCells, xPosition, yPosition, occupied, deathTime, totalNumberOfTimeSteps):

    cellNumber = 0

    # seed initial cells
    # space cells out evenly along the length of the capillary
    for x in range(int(0.5*xSteps/numberOfCells), xSteps, int(xSteps/numberOfCells)):
        xPosition[cellNumber][0] = x
        yPosition[cellNumber][0] = 0
        occupied[0][x] += 1

        # set the death time of all of the cells to the max time the simulation will run
        deathTime[cellNumber] = totalNumberOfTimeSteps - 1
        cellNumber += 1

    return
