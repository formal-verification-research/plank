# Description
# anastomosis kills the EC if they run into a different cell line


# Function
def anastomosis(anastomotic, y_position, workspace, file_deaths, death_time, occupied, cell, current_time_step,
                x_position, cell_line):

    # Set some position values for easier indexing
    y = y_position[cell][current_time_step]
    x = x_position[cell][current_time_step]
    y1 = y_position[cell][current_time_step + 1]
    x1 = x_position[cell][current_time_step + 1]

    workspace[y][x] = cell_line[cell]  # Set the EC to the current cell_line so it doesn't kill itself

    # anastomotic is turned on
    if anastomotic is True:

        # It must be outside the parent capillary
        if y1 > 0:

            # Keeps it from killing itself on its own cell_line or ECM
            if workspace[y1][x1] != 0 and workspace[y1][x1] != cell_line[cell]:

                # Kill it; it ran into a different capillary
                file_deaths.write("\n\nTime: " + str(current_time_step) + "\n")
                file_deaths.write("Cell: " + str(cell) + "\n")
                file_deaths.write("Cell ran into another capillary" + "\n")
                death_time[cell] = current_time_step + 1
                occupied[1][y1][x1] -= 1
                workspace[y1][x1] = cell_line[cell]

    # Set the EC back so it can be distinguished from the capillary
    workspace[y_position[cell][current_time_step + 1]][x_position[cell][current_time_step + 1]] = 100

    return workspace, file_deaths, death_time, occupied
