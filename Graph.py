# Description
'''
The purpose of Graph is to turn the difference substrate matrices and the workspace matrix into a graph for better
visual interpretation of the results
'''


# Imports
from numpy import linspace
from matplotlib import pyplot
from mpl_toolkits import mplot3d
from matplotlib import cm


# Function
def createGraph(ySubstrate, xSteps, vegf, pedf, fibronectin, protease, xVector, yVector, workspace, currentTimeStep):

    # Create the vegf height data for the 3D graph
    VEGFzVector = []
    for y in range(ySubstrate):
        if y % 2 == 0:
            for x in range(xSteps-1):
                VEGFzVector.append(vegf[y][x])
        else:
            for x in range(xSteps):
                VEGFzVector.append(vegf[y][x])

    # Create the pedf height data for the 3D graph
    PEDFzVector = []
    for y in range(ySubstrate):
        if y % 2 == 0:
            for x in range(xSteps - 1):
                PEDFzVector.append(pedf[y][x])
        else:
            for x in range(xSteps):
                PEDFzVector.append(pedf[y][x])

    # Create the vegf - pedf height data for the 3D graph
    VPzVector = []
    for y in range(ySubstrate):
        if y % 2 == 0:
            for x in range(xSteps-1):
                VPzVector.append(vegf[y][x] - pedf[y][x])
        else:
            for x in range(xSteps):
                VPzVector.append(vegf[y][x] - pedf[y][x])

    # Create the fibronectin height data for the 3D graph
    FIBRONECTINzVector = []
    for y in range(ySubstrate):
        if y % 2 == 0:
            for x in range(xSteps-1):
                FIBRONECTINzVector.append(fibronectin[y][x])
        else:
            for x in range(xSteps):
                FIBRONECTINzVector.append(fibronectin[y][x])

    # Create the protease height data for the 3D graph
    PROTEASEzVector = []
    for y in range(ySubstrate):
        if y % 2 == 0:
            for x in range(xSteps - 1):
                PROTEASEzVector.append(protease[y][x])
        else:
            for x in range(xSteps):
                PROTEASEzVector.append(protease[y][x])

    # Create the overall file
    fileName = "Time = " + str(currentTimeStep) + ".pdf"
    file = open(fileName, "a+")
    fig = pyplot.figure(figsize=pyplot.figaspect(4.0))

    # Create the EC color map
    ax = fig.add_subplot(7, 1, 1)
    ax.imshow(workspace, cmap='Purples')

    # Create the vegf 3D graph
    ax = fig.add_subplot(7, 1, 2, projection='3d')
    ax.plot_trisurf(xVector, yVector, VEGFzVector, cmap='hot', edgecolor='none')

    # Create the pedf 3D graph
    ax = fig.add_subplot(7, 1, 3, projection='3d')
    ax.plot_trisurf(xVector, yVector, PEDFzVector, cmap='cool', edgecolor='none')

    # Create the vp 3D graph
    ax = fig.add_subplot(7, 1, 4, projection='3d')
    ax.plot_trisurf(xVector, yVector, VPzVector, cmap='summer', edgecolor='none')

    # Create the fibronectin 3D graph
    ax = fig.add_subplot(7, 1, 5, projection='3d')
    ax.plot_trisurf(xVector, yVector, FIBRONECTINzVector, cmap='plasma', edgecolor='none')

    # Create the protease 3D graph
    ax = fig.add_subplot(7, 1, 6, projection='3d')
    ax.plot_trisurf(xVector, yVector, PROTEASEzVector, cmap='viridis', edgecolor='none')

    # 2D VEGF graph
    ax = fig.add_subplot(7, 1, 7)
    ax.imshow(vegf, cmap='hot')

    # Save and close the file so the program can run unattended
    pyplot.savefig(fileName)
    file.close()

    return
