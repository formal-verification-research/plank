"""
Description:
The simulation file takes the model created by the startup file and runs the simulation until completion.
"""


# Imports
from random import random
from prob_stay import prob_stay
from prob_move import prob_move
from move import move
from proliferation import proliferation
from anastomosis import anastomosis
from update_vegf import update_vegf
from update_pedf import update_pedf
from update_fib import update_fib
from update_pro import update_pro
from graph import graph
from parameter_vault import x_steps, threshold, graph_time


# Function
def simulation(x_length, y_steps, y_substrate, file_events, total_time, total_number_time_steps,
               k, h, lam, x_vector, y_vector, x_position, y_position, death_time, birth_time, divide_time, vegf, pedf,
               pro, fib, vegf_old, pedf_old, pro_old, fib_old, model, ec, ec_old, density_cap, density_ecm,
               cell_lineage, cell_tracker, child, number_of_cells):

    # Start the 'for-loop' that will take the simulation through the time steps
    for current_time_step in range(total_number_time_steps - 1):

        # Copy the current EC array into the previous one
        for x in range(x_steps):
            for y in range(y_steps):
                ec_old[y][x] = ec[y][x]

        # Start the 'for-loop' that will cycle through each EC for each time step
        for cell in range(number_of_cells):

            # Use the tip_cell function to run the processes that happen to a single EC for a single time step
            # At first assume the EC will not move
            x = x_position[cell][current_time_step]
            y = y_position[cell][current_time_step]
            x_position[cell][current_time_step+1] = x
            y_position[cell][current_time_step+1] = y

            # Perform the following actions if the EC is alive, which means the death time is still set at the end
            # Perform the following actions only if the EC is in the ECM, not in the parent blood vessel
            if death_time[cell] == total_number_time_steps - 1 and y > 0:

                # Add the current time step, Y and X EC positions to the cell tracking vector to find speeds, distances
                cell_tracker[cell][0].append(current_time_step)
                cell_tracker[cell][1].append(x)
                cell_tracker[cell][2].append(y)

                # Kill the EC if it reaches the RPE layer so it doesn't start backtracking
                if y == y_steps - 1:
                    death_time[cell] = current_time_step
                    ec[y][x] -= 1

                # Run the proliferation function to see if the EC will possibly divide or die
                number_of_cells = \
                    proliferation(y_steps, file_events, total_number_time_steps, k, x_position, y_position,
                                    death_time, birth_time, divide_time, pro, fib, pro_old, fib_old, model, ec,
                                  cell_lineage, number_of_cells, child, cell, current_time_step, x, y)

            # Find the probability that the EC will stay or move. Generate the random number that decides the EC outcome
            if death_time[cell] == total_number_time_steps - 1:
                stay = prob_stay(y, lam, k)
                left, T = prob_move(x, y, 0, pro, fib, vegf, y_steps, lam, k)
                right, T = prob_move(x, y, 1, pro, fib, vegf, y_steps, lam, k)
                up, T = prob_move(x, y, 2, pro, fib, vegf, y_steps, lam, k)
                random_num = random()

                # Perform the following action only if the EC is in the parent blood vessel
                if y == 0:

                    # Check the Fibronectin in the capillary wall to see if the EC can escape the parent blood vessel
                    if x == 0:
                        fib_cap = fib[0][0]
                    elif x == x_steps - 1:
                        fib_cap = fib[0][x_steps-2]
                    else:
                        fib_cap = (fib[0][x-1] + fib[0][x]) / 2
                    if fib_cap < threshold:
                        random_num = 2
                        file_events.write("\n\nTime: " + str(current_time_step) + "\n")
                        file_events.write("Cell: " + str(cell) + "\n")
                        file_events.write("Left the capillary" + "\n")

                # Use the move function to carry out the move decision made by the random number generator
                move(cell, current_time_step, stay, left, right, up, random_num, y_position, x_position, ec)

                # Use the anastomosis function to merge any capillaries if an EC happens to hit another capillary
                anastomosis(y_position, model, file_events, death_time, ec, cell, current_time_step, x_position,
                                  cell_lineage, total_number_time_steps)

        # As a shortcut you can end the simulation early if all of the EC had already died
        deaths = 0
        for cell in range(number_of_cells):
            if death_time[cell] != total_number_time_steps - 1:
                deaths += 1
        if deaths == number_of_cells:
            print("All of the EC have died")
            break

        # Update the substrate matrices, this is the part of the simulation that takes so long for the PC to run
        vegf = update_vegf(y_substrate, density_cap, density_ecm, ec_old, vegf, vegf_old, k, h, x_length)
        pedf = update_pedf(y_substrate, density_cap, density_ecm, ec_old, pedf, pedf_old, k, h, x_length)
        fib, fib_old = \
            update_fib(y_substrate, density_cap, ec_old, fib, fib_old, k, pro, h)
        pro, pro_old = \
            update_pro(y_substrate, density_cap, density_ecm, ec_old, pro, pro_old, k, vegf, pedf)

        # Print a notice to the screen so we know how far along the simulation is
        completion = current_time_step / total_number_time_steps * 100
        alive_cells = number_of_cells - deaths
        print("Time Progress = " + str(round(completion, 3)) + "%")
        print("Alive EC: " + str(alive_cells))
        print("Time step: " + str(current_time_step) + "\n")

        # Create a graph to show the progression of the simulation every specified amount of time
        if current_time_step % graph_time == 0 or current_time_step == 0:
            graph(y_substrate, vegf, pedf, fib, pro, x_vector, y_vector, model, current_time_step,
                  total_number_time_steps, total_time)

    graph(y_substrate, vegf, pedf, fib, pro, x_vector, y_vector, model, current_time_step, total_number_time_steps,
          total_time)

    return
