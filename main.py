"""
Description:
The main file starts the model code and runs all other files. It takes the input parameters, builds the model, runs
the simulation, and creates all graphs and other data outputs. This code predicts the movement of EC cells exposed to a
VEGF source, and then predicts changes caused by the introduction of PEDF to the model.
"""


# Imports
from startup import startup
from simulation import simulation


# Function
def main():

    # Create the model using the startup file
    start_time, x_length, y_length, y_steps, y_substrate, file_events, total_time, total_number_time_steps, k, h, \
    lam, x_vector, y_vector, x_position, y_position, death_time, birth_time, divide_time, vegf, pedf, pro, \
    fib, vegf_old, pedf_old, pro_old, fib_old, model, ec, ec_old, density_cap, density_ecm, cell_lineage, \
    cell_tracker, cell_number, threshold, child \
        = startup()




    # Run the simulation
    simulation(start_time, x_length, y_length, y_steps, y_substrate, file_events, total_time, total_number_time_steps,
               k, h, lam, x_vector, y_vector, x_position, y_position, death_time, birth_time, divide_time, vegf, pedf,
               pro, fib, vegf_old, pedf_old, pro_old, fib_old, model, ec, ec_old, density_cap, density_ecm,
               cell_lineage, cell_tracker, cell_number, threshold, child)

    # Finish up and the end
    end()

    return


# Call the main function
main()
