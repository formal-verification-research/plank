# Description
# anastomosis kills the EC if it runs into a different cell line


# Function
def anastomosis(anastomotic, y_position, workspace, file_deaths, death_time, occupied, cell, current_time_step,
                x_position, cell_line):

    # Set some position values for easier indexing within this function
    y = y_position[cell][current_time_step]
    x = x_position[cell][current_time_step]
    y1 = y_position[cell][current_time_step + 1]
    x1 = x_position[cell][current_time_step + 1]

    # Set the EC to the current cell_line so it doesn't kill itself
    workspace[y][x] = cell_line[cell]

    if anastomotic is True:  # anastomotic is turned on
        if y1 > 0:  # It must be outside the parent capillary
            # Keeps it from killing itself on its own cell_line or ECM
            if workspace[y1][x1] != 0 and workspace[y1][x1] != cell_line[cell]:
                # Kill it; it ran into a different capillary
                file_deaths.write("\n\nTime: " + str(current_time_step) + "\n")
                file_deaths.write("Cell: " + str(cell) + "\n")
                file_deaths.write("Cell ran into another capillary" + "\n")
                death_time[cell] = current_time_step + 1
                occupied[y1][x1] -= 1
                workspace[y1][x1] = cell_line[cell]

    # Set the EC back so it can be distinguished from the capillary
    workspace[y1][x1] = 100

    return
