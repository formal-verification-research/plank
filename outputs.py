"""
Description:
The outputs file records the values of directness, distance, and speed to be compared to in-vitro experiments
"""


# Imports
from numpy import sqrt
from parameter_vault import x_steps, L, time_step_duration


# Function
def outputs(cell_tracker, x_length):

    # Create the file and get the individual EC data only if they traveled more than 1 space
    file_outputs = open("Outputs.txt", "w")
    for cell in range(len(cell_tracker)):
        if len(cell_tracker[cell][0]) < 2:
            continue

        # Create the distance vector and find distances with the pythagorean formula
        distances = []
        for time in range(1, len(cell_tracker[cell][0])):
            distances.append(sqrt((cell_tracker[cell][1][time] - cell_tracker[cell][1][time - 1]) ** 2
                                  + (cell_tracker[cell][2][time] - cell_tracker[cell][2][time - 1]) ** 2))

        # Calculate all the different output information that we are interested in
        distance_accumulated = sum(distances)
        distance_euclidean = sqrt((cell_tracker[cell][1][-1] - cell_tracker[cell][1][0]) ** 2
                                  + (cell_tracker[cell][2][-1] - cell_tracker[cell][2][0]) ** 2)
        directness = distance_euclidean / distance_accumulated
        FMI_y = (cell_tracker[cell][2][-1] - cell_tracker[cell][2][0]) / distance_accumulated
        FMI_x = (cell_tracker[cell][1][-1] - cell_tracker[cell][1][0]) / distance_accumulated
        dist_accum_dimensionalized = distance_accumulated * (x_length / (x_steps-1)) * L * 1000  # um
        dist_euclid_dimensionalized = distance_euclidean * (x_length / (x_steps-1)) * L * 1000  # um
        time_dimensionalized = (cell_tracker[cell][0][-1] - cell_tracker[cell][0][0]) * time_step_duration / 60  # min
        speed_euclidean = dist_euclid_dimensionalized / time_dimensionalized
        speed_accumulated = dist_accum_dimensionalized / time_dimensionalized

        # Add the information outputs to the file
        file_outputs.write('Cell Number: ' + str(cell) + "\r\n")
        file_outputs.write('Directness: ' + str(directness) + "\r\n")
        file_outputs.write('Accumulated Distance (um): ' + str(dist_accum_dimensionalized) + "\r\n")
        file_outputs.write('Accumulated Speed (um/min): ' + str(speed_accumulated) + "\r\n")
        file_outputs.write('Euclidean Distance (um): ' + str(dist_euclid_dimensionalized) + "\r\n")
        file_outputs.write('Euclidean Speed (um/min): ' + str(speed_euclidean) + "\r\n")
        file_outputs.write('Time (min): ' + str(time_dimensionalized) + "\r\n\r\n")

    # Make sure to close the file so all the information will save
    file_outputs.close()

    return
