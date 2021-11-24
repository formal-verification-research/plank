# Description
# data_outputs records the output values directness, distance, and velocity to be compared to in-vitro experiments


# Imports
from numpy import sqrt


# Function
def data_outputs(cell_tracker, x_length, x_steps, L):

    file_outputs = open("Data_Outputs.txt", "w")  # Create the file
    for cell in range(len(cell_tracker)):  # Get individual data for each active EC
        if len(cell_tracker[cell][0]) < 2:  # Skip if EC didn't go more than 1 space
            continue

        distances = []  # Initialize vector to collect distance moved each time step
        for time in range(1, len(cell_tracker[cell][0])):  # Calculate distance with pythagorean formula
            distances.append(sqrt((cell_tracker[cell][1][time] - cell_tracker[cell][1][time - 1]) ** 2
                                  + (cell_tracker[cell][2][time] - cell_tracker[cell][2][time - 1]) ** 2))

        # Calculate outputs
        distance_accumulated = sum(distances)
        distance_euclidean = sqrt((cell_tracker[cell][1][-1] - cell_tracker[cell][1][0]) ** 2
                                  + (cell_tracker[cell][2][-1] - cell_tracker[cell][2][0]) ** 2)
        directness = distance_euclidean / distance_accumulated
        FMI_y = (cell_tracker[cell][2][-1] - cell_tracker[cell][2][0]) / distance_accumulated
        FMI_x = (cell_tracker[cell][1][-1] - cell_tracker[cell][1][0]) / distance_accumulated
        dist_accum_dimensionalized = distance_accumulated * (x_length / (x_steps-1)) * L * 1000  # um
        dist_euclid_dimensionalized = distance_euclidean * (x_length / (x_steps-1)) * L * 1000  # um
        time_dimensionalized = (cell_tracker[cell][0][-1] - cell_tracker[cell][0][0]) * 8 / 60  # min, 8s = time step
        velocity_euclidean = dist_euclid_dimensionalized / time_dimensionalized
        velocity_accumulated = dist_accum_dimensionalized / time_dimensionalized

        # Add outputs to the file
        file_outputs.write('Cell Number: ' + str(cell) + "\r\n")
        file_outputs.write('Directness: ' + str(directness) + "\r\n")
        file_outputs.write('Accumulated Distance (um): ' + str(dist_accum_dimensionalized) + "\r\n")
        file_outputs.write('Accumulated Velocity (um/min): ' + str(velocity_accumulated) + "\r\n")
        file_outputs.write('Euclidian Distance (um): ' + str(dist_euclid_dimensionalized) + "\r\n")
        file_outputs.write('Euclidian Velocity (um/min): ' + str(velocity_euclidean) + "\r\n")
        file_outputs.write('Time (min): ' + str(time_dimensionalized) + "\r\n\r\n")

    file_outputs.close()  # Close the file

    return
