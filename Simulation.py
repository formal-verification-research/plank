# Description
'''
The purpose of Simulation is to run the movement of the EC from the starting point until they all die
'''


# Imports
# Library Imports
from random import random
import numpy as np
import math
from scipy import spatial
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib import pyplot
import matplotlib.backends.backend_pdf

# File Imports
from Heaviside import heaviside
from pStay import pStay
from pMove import pMove
from Move import move
from UpdateProtease import updateProtease
from UpdateVEGF import updateVEGF
from UpdatePEDF import updatePEDF
from UpdateFibronectin import updateFibronectin
from Graph import createGraph
from proliferation import proliferation
from anastomosis import anastomosis

# Values Imports
from Values import K18, K20, K25, K26, m1

# Function
def simulation(totalNumberOfTimeSteps, xSteps, ySteps, occupied, occupiedOld, numberOfCells, xPosition, yPosition,
               deathTime, protease, proteaseOld, densityScale, lam, k, fibronectin, vegf, ySubstrate, vegfOld,
               tolerance, h, xLength, fibronectinOld, xVector, yVector, maxCellsAllowed, birthTime, divideTime,
               fibronectinThreshold, graphTime, child, anastomotic, cellLine, fileDeaths, fileDivisions, cellTrackingVector,
               pedf, pedfOld, injection):

    # Create the EC matrix
    workspace = np.zeros((ySteps, xSteps))

    # Cycle through time steps
    for currentTimeStep in range(totalNumberOfTimeSteps - 1):

        # Copy occupied into occupiedOld
        for x in range(xSteps):
            for y in range(ySteps):
                occupiedOld[y][x] = occupied[y][x]

        # Cycle through cells
        for cell in range(numberOfCells):
            # At first assume no movement
            x = xPosition[cell][currentTimeStep]
            y = yPosition[cell][currentTimeStep]
            xPosition[cell][currentTimeStep + 1] = x
            yPosition[cell][currentTimeStep + 1] = y

            # If cell has left the capillary
            # Determine if the cell has divided or died
            if deathTime[cell] == totalNumberOfTimeSteps - 1 and y > 0:

                # Add the time, x, y to the cell tracking vector
                cellTrackingVector[cell][0].append(currentTimeStep)
                cellTrackingVector[cell][1].append(x)
                cellTrackingVector[cell][2].append(y)

                # The cell dies/leaves simulation if it reaches the tumour
                if y == ySteps - 1:
                    deathTime[cell] = currentTimeStep
                    occupied[y][x] -= 1

                # Calculate average protease values at currentTimeStep and currentTimeStep-1
                # Surrounding mesh points usually equals 4 unless it is a boundary condition
                proteaseMinus0 = 0
                proteaseMinus1 = 0
                surroundingPoints = 0

                # Left
                if x > 0:
                    proteaseMinus0 = proteaseMinus0 + protease[2 * y][x - 1]
                    proteaseMinus1 = proteaseMinus1 + proteaseOld[2 * y][x - 1]
                    surroundingPoints += 1

                # Right
                if x < xSteps - 1:
                    proteaseMinus0 = proteaseMinus0 + protease[2 * y][x]
                    proteaseMinus1 = proteaseMinus1 + proteaseOld[2 * y][x]
                    surroundingPoints += 1

                # Up
                if y > 0:
                    proteaseMinus0 = proteaseMinus0 + protease[2 * y - 1][x]
                    proteaseMinus1 = proteaseMinus1 + proteaseOld[2 * y - 1][x]
                    surroundingPoints += 1

                # Down
                if y < ySteps - 1:
                    proteaseMinus0 = proteaseMinus0 + protease[2 * y + 1][x - 1]
                    proteaseMinus1 = proteaseMinus1 + proteaseOld[2 * y + 1][x - 1]
                    surroundingPoints += 1

                proteaseMinus0 = proteaseMinus0 / surroundingPoints
                proteaseMinus1 = proteaseMinus1 / surroundingPoints

                # Find G, used for divide probabilities
                G = K25 * (math.exp(-K26 * proteaseMinus0 ** m1) * (1 - K26 * m1 * proteaseMinus0 ** m1)) / \
                    (1 + K25 * proteaseMinus0 * math.exp(-K26 * proteaseMinus0 ** m1))

                # Find the protease dependant term for divide and death probabilities
                proteaseDependent = G * (proteaseMinus0 - proteaseMinus1) / k

                # Find divide and death probabilities
                if proteaseDependent >= 0:
                    divideProbability = (k * K18 + G * (proteaseMinus0 - proteaseMinus1)) * heaviside(currentTimeStep - divideTime[cell] - child * totalNumberOfTimeSteps)
                    deathProbability = k * K20
                else:
                    divideProbability = k * K18 * heaviside(currentTimeStep - divideTime[cell] - child * totalNumberOfTimeSteps)
                    deathProbability = k * K20 - G * (proteaseMinus0 - proteaseMinus1)

                proliferation(deathProbability, deathTime, cell, currentTimeStep, occupied, y, x, workspace, cellLine,
                              fileDeaths, divideProbability, numberOfCells, maxCellsAllowed, xPosition, yPosition,
                              totalNumberOfTimeSteps, birthTime, divideTime, fileDivisions)

            # Determine if the cell moves and where
            if deathTime[cell] == totalNumberOfTimeSteps - 1:
                stay = pStay(y, lam, k)
                left, T = pMove(x, y, 0, protease, fibronectin, vegf, pedf, xSteps, ySteps, lam, k)
                right, T = pMove(x, y, 1, protease, fibronectin, vegf, pedf, xSteps, ySteps, lam, k)
                up, T = pMove(x, y, 2, protease, fibronectin, vegf, pedf, xSteps, ySteps, lam, k)
                randomNumber = random()

                # Check if cell can escape the capillary
                if y == 0:
                    if x == 0:
                        fibronectinCapillary = fibronectin[0][0]
                    elif x == xSteps - 1:
                        fibronectinCapillary = fibronectin[0][xSteps - 2]
                    else:
                        fibronectinCapillary = (fibronectin[0][x - 1] + fibronectin[0][x]) / 2
                    if fibronectinCapillary < fibronectinThreshold:
                        randomNumber = 2

                # Move the EC
                move(cell, currentTimeStep, stay, left, right, up, randomNumber, yPosition, xPosition, occupied)

                # Anastomosis and workspace
                anastomosis(anastomotic, yPosition, workspace, fileDeaths, deathTime, occupied, cell, currentTimeStep,
                    xPosition, cellLine)

        # Find out when all the EC have died, and end the program early
        deaths = 0
        for cell in range(numberOfCells):
            if deathTime[cell] != totalNumberOfTimeSteps - 1:
                deaths += 1
        if deaths == numberOfCells:
            print("All of the EC have died")
            break

        updateVEGF(ySubstrate, xSteps, densityScale, occupiedOld, vegf, vegfOld, k, tolerance, h, xLength)
        #if currentTimeStep >= injection:
            #updatePEDF(ySubstrate, xSteps, densityScale, occupiedOld, pedf, pedfOld, k, tolerance, h, xLength)
        updateFibronectin(ySubstrate, xSteps, densityScale, occupiedOld, fibronectin, fibronectinOld, k, protease,
                          tolerance, h)
        updateProtease(ySubstrate, xSteps, densityScale, occupiedOld, protease, proteaseOld, k, vegfOld, pedfOld)

        print("Current Time Step = " + str(currentTimeStep))

        if currentTimeStep % graphTime == 0:
            createGraph(ySubstrate, xSteps, vegf, pedf, fibronectin, protease, xVector, yVector, workspace, currentTimeStep)

    createGraph(ySubstrate, xSteps, vegf, pedf, fibronectin, protease, xVector, yVector, workspace, currentTimeStep)

    return
