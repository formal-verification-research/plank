# Description
# end closes the files and prints any remaining values needed


# Imports
from data_outputs import data_outputs
import time

# Function
def end(cell_tracker, x_length, x_steps, L, start_time, file_events):

    # Get the data outputs
    data_outputs(cell_tracker, x_length, x_steps, L)

    print("Time it took to run:\n" + str((time.time()-start_time)/60/60) + "hours")  # How much time it took to run
    file_events.close()  # Close the files

    return
