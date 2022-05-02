"""
Description:
The anastomosis file merges two capillary trails into a loop if a EC happens to collide with a different capillary
trail. It does this by 'killing' = (changing the death time to the current time step)' of the EC that entered the
capillary of another. The file judges this based on the model array, which holds the location of all the EC,
capillaries, and ECM. Currently, this file will only kill an EC if it enters the capillary of a different cell lineage.
If the EC runs into a capillary of the same cell lineage or another EC, it will allow both EC to continue.
"""


# Imports
from parameter_vault import anastomotic


# Function
def anastomosis(y_position, model, file_events, death_time, ec, cell, current_time_step, x_position, cell_lineage,
                total_number_time_steps):

    # Set the lengthy index value to a simple variable to make indexing easier within this function
    y = y_position[cell][current_time_step]
    x = x_position[cell][current_time_step]
    y1 = y_position[cell][current_time_step+1]
    x1 = x_position[cell][current_time_step+1]

    # Set the previous point where the EC was to be a capillary trail mark on the model now
    model[y][x] = cell_lineage[cell]

    # The anastomosis part is set to run only if the anastomotic toggle is turned on and the EC is in the ECM
    if anastomotic is True:
        if y1 > 0:
            if model[y1][x1] != 0 and model[y1][x1] != cell_lineage[cell] and model[y1][x1] != 100:
                file_events.write("\n\nTime: " + str(current_time_step) + "\n")
                file_events.write("Cell: " + str(cell) + "\n")
                file_events.write("Cell ran into another capillary" + "\n")
                death_time[cell] = current_time_step + 1
                ec[y1][x1] -= 1
                model[y1][x1] = cell_lineage[cell]

    # Set the EC back so it can be distinguished from the capillary on the model array
    if death_time[cell] == total_number_time_steps - 1:
        model[y1][x1] = 100

    return
