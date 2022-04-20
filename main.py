# Description
# main starts the program and runs all other files. This experiment predicts the movement of EC cells exposed to a
# vegf source, and then predicts changes caused by the introduction of pedf to the model.


# Imports
from start import terms
from start import arrays
from start import parent_vessel
from simulation import simulation
from end import end
from parameter_vault import x_steps, number_of_cells, tolerance, max_cells_allowed, graph_time, total_time, \
    total_number_time_steps, threshold, child, anastomotic, L


# Function
def main():

    # Create the starting terms and parameters
    start_time, x_length, y_length, y_steps, y_substrate, file_events, k, h, lam, x_vector, y_vector \
        = terms(x_steps, total_time, total_number_time_steps)

    # Create the starting arrays
    x_position, y_position, death_time, birth_time, divide_time, vegf, pedf, pro, fib, vegf_old, pedf_old, pro_old, \
        fib_old, workspace \
        = arrays(max_cells_allowed, total_number_time_steps, y_substrate, x_steps, y_steps)

    # Create the starting capillary
    occupied, occupied_old, density_scale, cell_line, cell_tracker, cell_number \
        = parent_vessel(y_steps, x_steps, number_of_cells, max_cells_allowed, x_position, y_position, death_time,
                  total_number_time_steps)

    # Run simulation
    simulation(total_number_time_steps, x_steps, y_steps, occupied, number_of_cells, x_position, y_position,
               death_time, pro, density_scale, lam, k, fib, vegf, y_substrate, tolerance, h, x_length, x_vector,
               y_vector, max_cells_allowed, birth_time, divide_time, threshold, graph_time, child, anastomotic,
               cell_line, file_events, cell_tracker, pedf, vegf_old, pedf_old, pro_old, fib_old, occupied_old,
               workspace, total_time)

    # Finish up and the end
    end(cell_tracker, x_length, x_steps, L, start_time, file_events)

    return


main()
