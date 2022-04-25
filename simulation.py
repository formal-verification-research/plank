"""
Description:
The simulation file takes the model created by the startup file and runs the simulation until completion.
"""


# Imports
from tip_cell import tip_cell
from update_vegf import update_vegf
from update_pedf import update_pedf
from update_fib import update_fib
from update_pro import update_pro
from graph import graph
from parameter_vault import x_steps


# Function
def simulation(x_length, y_steps, y_substrate, file_events, total_time, total_number_time_steps,
               k, h, lam, x_vector, y_vector, x_position, y_position, death_time, birth_time, divide_time, vegf, pedf,
               pro, fib, vegf_old, pedf_old, pro_old, fib_old, model, ec, ec_old, density_cap, density_ecm,
               cell_lineage, cell_tracker, threshold, child, graphing, number_of_cells):

    # Start the 'for-loop' that will take the simulation through the time steps
    for current_time_step in range(total_number_time_steps - 1):

        # Copy the current EC array into the previous one
        for x in range(x_steps):
            for y in range(y_steps):
                ec_old[y][x] = ec[y][x]

        # Start the 'for-loop' that will cycle through each EC for each time step
        for cell in range(number_of_cells):

            # Use the tip_cell function to run the processes that happen to a single EC for a single time step
            ec, ec_old, cell_tracker, death_time, model, file_events, number_of_cells, x_position, y_position, \
            cell_lineage, birth_time, divide_time \
                = tip_cell(y_steps, file_events, total_number_time_steps, k, lam, x_position, y_position, death_time,
                     birth_time, divide_time, vegf, pro, fib, pro_old, fib_old, model, ec, ec_old, cell_lineage,
                     cell_tracker, number_of_cells, threshold, child, cell, current_time_step)

        # As a shortcut you can end the simulation early if all of the EC had already died
        deaths = 0
        for cell in range(number_of_cells):
            if death_time[cell] != total_number_time_steps - 1:
                deaths += 1
        if deaths == number_of_cells:
            print("All of the EC have died")
            break

        # Update the substrate matrices, this is the part of the simulation that takes so long for the PC to run
        vegf, vegf_old = \
            update_vegf(y_substrate, density_cap, density_ecm, ec_old, vegf, vegf_old, k, h, x_length)
        pedf, pedf_old = \
            update_pedf(y_substrate, density_cap, density_ecm, ec_old, pedf, pedf_old, k, h, x_length)
        fib, fib_old = \
            update_fib(y_substrate, density_cap, ec_old, fib, fib_old, k, pro, h)
        pro, pro_old = \
            update_pro(y_substrate, density_cap, density_ecm, ec_old, pro, pro_old, k, vegf_old, pedf_old)

        # Print a notice to the screen so we know how far along the simulation is
        completion = current_time_step / total_number_time_steps * 100
        alive_cells = number_of_cells - deaths
        print("Time Progress = " + str(round(completion, 3)) + "%")
        print("Alive EC: " + str(alive_cells))
        print("Time step: " + str(current_time_step) + "\n\n")

        # Create a graph to show the progression of the simulation every specified amount of time
        if current_time_step % graphing == 0 or current_time_step== 0:
            graph(y_substrate, vegf, pedf, fib, pro, x_vector, y_vector, model, current_time_step,
                  total_number_time_steps, total_time)

    graph(y_substrate, vegf, pedf, fib, pro, x_vector, y_vector, model, current_time_step, total_number_time_steps,
          total_time)

    return cell_tracker, file_events
