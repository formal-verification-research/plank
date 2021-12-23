# Description
# main starts the program and runs all other files. This experiment predicts the movement of EC cells exposed to a
# vegf source, and then predicts changes caused by the introduction of pedf to the model.


# Imports
from start import terms
from start import arrays
from start import parent_vessel
from simulation import simulation
from end import end


# This control_center lists the variables that are allowed to be modified by the user. The variables are passed into
# the main function. These variables are changeable because the user might want to test different ideas.
from the_vault import L  # L can be changed; L must be in the_vault
x_steps = 21  # Nodes in the x domain
number_of_cells = 5  # How many EC to start with
tolerance = 0.001  # Accuracy tolerance of the substrate updaters
max_cells_allowed = 100  # How many total cells are allowed in the experiment
graph_time = 200  # How often a graph is created in amount of time steps
total_time = 0.06912  # Dimension-less time the simulation lasts. 48 hours = 0.06912, plank pg 150.
total_number_time_steps = 21600  # How many time steps
threshold = 0.6  # Level that fib has to drop to before EC can leave the parent capillary
child = 450  # How long between divisions in time steps
anastomotic = True  # Anastomosis


# Function
def main(x_steps, number_of_cells, tolerance, max_cells_allowed, graph_time, total_time, total_number_time_steps,
         threshold, child, anastomotic, L):

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


main(x_steps, number_of_cells, tolerance, max_cells_allowed, graph_time, total_time, total_number_time_steps,
     threshold, child, anastomotic, L)
