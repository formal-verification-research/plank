"""
Description:
The close file closes the open files so they can have the data included. It also prints any remaining values needed
"""


# Imports
from outputs import outputs
import time


# Function
def close(cell_tracker, x_length, start_time, file_events):

    # Run the data_outputs function to collect info from the simulation such as speed and distance
    outputs(cell_tracker, x_length)

    # Print how much time it took the computer to run the code and close the events file
    print("Time it took to run:\n" + str((time.time()-start_time)/60/60) + "hours")
    file_events.close()

    return
