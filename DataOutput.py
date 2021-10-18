# Description
'''
The purpose of Data Output is to record the output values to be compared to in vitro experiments
'''


# Imports
import numpy as np


# Function
# File takes in the cell positions at each time (x and y)
# Calculates the FMIx, FMIy, Directness, Euclidean distance, Velocity
def dataOutput(cellTrackingVector, xLength, xSteps, totalTime, totalNumberOfTimeSteps):

    for cell in range(len(cellTrackingVector)):
        # Don't measure empty vectors
        if len(cellTrackingVector[cell][0]) == 0:
            continue

        # Loop through data and remove the data when the cells are still in the parent vessel
        for i in range(len(cellTrackingVector[cell][0])):
            if cellTrackingVector[cell][2][i] != 0:
                index = i
                break
        del cellTrackingVector[cell][0][:index]
        del cellTrackingVector[cell][1][:index]
        del cellTrackingVector[cell][2][:index]

        distances = []    # Initialize vector to collect distance moved each time step
        for time in range(1, len(cellTrackingVector[cell][0])):
            # Calculate the distance at each time step using the pythagorean formula
            distances.append(np.sqrt((cellTrackingVector[cell][1][time] - cellTrackingVector[cell][1][time-1])**2 + \
                (cellTrackingVector[cell][2][time] - cellTrackingVector[cell][2][time-1])**2))

        # Calculate outputs
        dist_accum = sum(distances)
        dist_euclid = np.sqrt((cellTrackingVector[cell][1][-1] - cellTrackingVector[cell][1][0])**2 + \
            (cellTrackingVector[cell][2][-1] - cellTrackingVector[cell][2][0])**2)
        directness = dist_euclid / dist_accum
        FMI_y = (cellTrackingVector[cell][2][-1] - cellTrackingVector[cell][2][0]) / dist_accum
        FMI_x = (cellTrackingVector[cell][1][-1] - cellTrackingVector[cell][1][0]) / dist_accum

        dist_accum_nondimensionalized = dist_accum * (xLength/xSteps) * 1000 # Gives distance in um
        dist_euclid_nondimensionalized = dist_euclid * (xLength/xSteps) * 1000 # Gives distance in um

        time_nondimensionalized = (cellTrackingVector[cell][0][-1] - cellTrackingVector[cell][0][0]) * 8 / 60 # each time step is 8 seconds
        velocity_euclid_nondimensionalized = dist_euclid_nondimensionalized / time_nondimensionalized
        velocity_accum_nondimensionalized = dist_accum_nondimensionalized / time_nondimensionalized

        # Create a file
        fileOutputs = open("Data_Outputs{}.txt".format(cell), "w")
        fileOutputs.write('Cell Number: ' + str(cell) + "\r\n")
        fileOutputs.write('Accumulated Distance: ' + str(dist_accum) + "\r\n")
        fileOutputs.write('Euclidian Distance: ' + str(dist_euclid) + "\r\n")
        fileOutputs.write('Directness: ' + str(directness) + "\r\n")
        fileOutputs.write('FMI_y: ' + str(FMI_y) + "\r\n")
        fileOutputs.write('FMI_x: ' + str(FMI_x) + "\r\n")
        fileOutputs.write('Accumulated Distance (um): ' + str(dist_accum_nondimensionalized) + "\r\n")
        fileOutputs.write('Euclidian Distance (um): ' + str(dist_euclid_nondimensionalized) + "\r\n")
        fileOutputs.write('Euclidian Velocity (um/min): ' + str(velocity_euclid_nondimensionalized) + "\r\n")
        fileOutputs.write('Accumulated Velocity (um/min): ' + str(velocity_accum_nondimensionalized) + "\r\n")
        fileOutputs.write('Time (min): ' + str(time_nondimensionalized) + "\r\n")
        fileOutputs.close()

    return
