# General Description
'''
This program predicts the movement of EC cells exposed to a VEGF source, and then predicts changes caused by the
introduction of PEDF
'''


# Imports
# Library Imports
from numpy import zeros
from numpy import ones
import time
from numpy import linspace
from matplotlib import pyplot
from mpl_toolkits import mplot3d

# File Imports
from StartUp import startUp
from Lamda import lamda
from Simulation import simulation
from DataOutput import dataOutput


# Input variables, meaning variables the user is allowed to change based on requirements
# Size of the matrix
xSteps = 201
# How many EC to start with
numberOfCells = 5    # numCells must be less than or equal to xStep so there will be a place to put all the cells
# Accuracy tolerance of the substrate updaters
tolerance = 0.001
# How many total cells are allowed in the experiment
maxCellsAllowed = 100
# How often you want a graph created in amount of time steps
graphTime = 200
# How long do you want the simulation to last
# 694.4 hours = 28.9 days is the calculated time for which equals 1 when the simulation goes dimensionless
# So to get 48 hours, the total time must be 0.06912
totalTime = 0.06912
# How many time steps
totalNumberOfTimeSteps = 21600
# Level that fibronectin has to drop to before EC can leave the capillary
fibronectinThreshold = 0.6
# How long between divisions
child = 0.125

# Anastomosis
anastomotic = False

# Function
def main(xSteps, numberOfCells, tolerance, maxCellsAllowed, graphTime, totalTime, totalNumberOfTimeSteps,
         fibronectinThreshold, child, anastomotic):

    # Lock in the start time of the program
    startTime = time.time()

    # Dimensionless size of the matrix
    xLength = 1
    yLength = 0.5
    ySteps = int(xSteps * (yLength/xLength) + 0.5)    # need to add the 0.5 so that it rounds to the correct number
    ySubstrate = ySteps * 2 - 1

    # Create the death and division files
    fileDeaths = open("Tracking_Deaths.txt", "w")
    fileDivisions = open("Tracking_Divisions.txt", "w")

    # Density of cells in the capillary
    densityScale = xSteps / numberOfCells

    # Set up x position and y positions
    xPosition = zeros((maxCellsAllowed, totalNumberOfTimeSteps), dtype=int)
    yPosition = zeros((maxCellsAllowed, totalNumberOfTimeSteps), dtype=int)

    # Death, birth, and divide times
    deathTime = zeros(maxCellsAllowed, dtype=int)
    birthTime = zeros(maxCellsAllowed, dtype=int)
    divideTime = zeros(maxCellsAllowed, dtype=int)

    # Substrate matrices
    vegf = zeros((ySubstrate, xSteps))
    protease = zeros((ySubstrate, xSteps))
    fibronectin = ones((ySubstrate, xSteps))
    vegfOld = zeros((ySubstrate, xSteps))
    proteaseOld = zeros((ySubstrate, xSteps))
    fibronectinOld = ones((ySubstrate, xSteps))

    # To tell if that space has an EC in it or not
    occupied = zeros((ySteps, xSteps), dtype=int)
    occupiedOld = zeros((ySteps, xSteps), dtype=int)

    # Try and create different cell line colors
    cellLine = zeros(maxCellsAllowed)
    cl = 20
    for i in range(numberOfCells):
        cellLine[i] = cl
        cl += 15

    # Create k, the time step interval variable
    k = totalTime / totalNumberOfTimeSteps

    # Create h, the distance between mesh points variable
    h = xLength / xSteps

    # Create lam, the ratio between length and step size
    lam = lamda(xLength, h)

    # Create the cell tracking vector
    cellTrackingVector = [[[], [], []] for i in range(maxCellsAllowed)]

    # Create the x and y vectors used to create the 3D graphs
    xVector = []
    yVector = []

    # Place the values inside the vectors based on if they are an odd or even row
    xEvenVector = []
    for x in linspace(0.5, xSteps - 1.5, num=xSteps - 1):
        xEvenVector.append(x)

    xOddVector = []
    for x in linspace(0, xSteps - 1, num=xSteps):
        xOddVector.append(x)

    for y in range(ySubstrate):
        if y % 2 == 0:
            for i in xEvenVector:
                xVector.append(i)
                yVector.append(y)
        else:
            for j in xOddVector:
                xVector.append(j)
                yVector.append(y)

    # Run StartUp, which will place all the EC and set up the matrices
    startUp(xSteps, numberOfCells, xPosition, yPosition, occupied, deathTime, totalNumberOfTimeSteps)

    # Run Simulation, which will go until all cells are dead
    simulation(totalNumberOfTimeSteps, xSteps, ySteps, occupied, occupiedOld, numberOfCells, xPosition, yPosition,
               deathTime, protease, proteaseOld, densityScale, lam, k, fibronectin, vegf, ySubstrate, vegfOld,
               tolerance, h, xLength, fibronectinOld, xVector, yVector, maxCellsAllowed, birthTime, divideTime,
               fibronectinThreshold, graphTime, child, anastomotic, cellLine, fileDeaths, fileDivisions, cellTrackingVector)

    dataOutput(cellTrackingVector, xLength, xSteps, totalTime, totalNumberOfTimeSteps)

    # Print how much time it took to run
    print("Time it took to run:\n")
    print(str((time.time()-startTime)/60/60) + "hours")

    # Close the files
    fileDeaths.close()
    fileDivisions.close()

    return

main(xSteps, numberOfCells, tolerance, maxCellsAllowed, graphTime, totalTime, totalNumberOfTimeSteps,
     fibronectinThreshold, child, anastomotic)
