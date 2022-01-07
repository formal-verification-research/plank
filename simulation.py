# Description
# simulation runs the movement of the EC from the starting point until they all die and updates each time step


# Imports
from random import random
from prob_stay import prob_stay
from prob_move import prob_move
from move import move
from update_pro import update_pro
from update_vegf import update_vegf
from update_pedf import update_pedf
from update_fib import update_fib
from graph import graph
from proliferation import proliferation
from anastomosis import anastomosis


# Function
def simulation(total_number_time_steps, x_steps, y_steps, occupied, number_of_cells, x_position, y_position,
               death_time, pro, density_scale, lam, k, fib, vegf, y_substrate, tolerance, h, x_length, x_vector,
               y_vector, max_cells_allowed, birth_time, divide_time, threshold, graph_time, child, anastomotic,
               cell_line, file_events, cell_tracker, pedf, vegf_old, pedf_old, pro_old, fib_old, occupied_old, 
               workspace, total_time):

    # Cycle through time steps
    for current_time_step in range(total_number_time_steps - 1):

        # Copy occupied into occupied_old
        for x in range(x_steps):
            for y in range(y_steps):
                occupied_old[y][x] = occupied[y][x]

        # Cycle through cells
        for cell in range(number_of_cells):

            # At first assume no movement
            x = x_position[cell][current_time_step]
            y = y_position[cell][current_time_step]
            x_position[cell][current_time_step + 1] = x
            y_position[cell][current_time_step + 1] = y

            # If cell has left the capillary
            if death_time[cell] == total_number_time_steps - 1 and y > 0:

                # Add the time, x, y to the cell tracking vector
                cell_tracker[cell][0].append(current_time_step)
                cell_tracker[cell][1].append(x)
                cell_tracker[cell][2].append(y)

                # The cell dies/leaves simulation if it reaches the tumour
                if y == y_steps - 1:
                    death_time[cell] = current_time_step
                    occupied[y][x] -= 1

                number_of_cells = \
                    proliferation(death_time, cell, current_time_step, occupied, y, x, workspace, cell_line,
                                  file_events, number_of_cells, max_cells_allowed, x_position, y_position,
                                  total_number_time_steps, birth_time, divide_time, pro, pedf, fib, pro_old, pedf_old,
                                  fib_old, x_steps, y_steps, k, child)

            # Determine if the cell moves and where
            if death_time[cell] == total_number_time_steps - 1:
                stay = prob_stay(y, lam, k)
                left, T = prob_move(x, y, 0, pro, fib, vegf, x_steps, y_steps, lam, k)
                right, T = prob_move(x, y, 1, pro, fib, vegf, x_steps, y_steps, lam, k)
                up, T = prob_move(x, y, 2, pro, fib, vegf, x_steps, y_steps, lam, k)
                random_num = random()

                # Check if cell can escape the capillary
                if y == 0:
                    if x == 0:
                        fib_cap = fib[0][0]
                    elif x == x_steps - 1:
                        fib_cap = fib[0][x_steps - 2]
                    else:
                        fib_cap = (fib[0][x - 1] + fib[0][x]) / 2
                    if fib_cap < threshold:
                        random_num = 2

                # Move the EC
                move(cell, current_time_step, stay, left, right, up, random_num, y_position, x_position, occupied)

                # Anastomosis and workspace
                anastomosis(anastomotic, y_position, workspace, file_events, death_time, occupied, cell,
                            current_time_step, x_position, cell_line)

        # Find out when all the EC have died, and end the program early
        deaths = 0
        for cell in range(number_of_cells):
            if death_time[cell] != total_number_time_steps - 1:
                deaths += 1
        if deaths == number_of_cells:
            print("All of the EC have died")
            break

        update_vegf(y_substrate, x_steps, density_scale, occupied_old, vegf, vegf_old, k, tolerance, h, x_length)
        update_pedf(y_substrate, x_steps, density_scale, occupied_old, pedf, pedf_old, k, tolerance, h, x_length)
        update_fib(y_substrate, x_steps, density_scale, occupied_old, fib, fib_old, k, pro, tolerance, h)
        update_pro(y_substrate, x_steps, density_scale, occupied_old, pro, pro_old, k, vegf_old, pedf_old)

        print("Current Time Step = " + str(current_time_step))

        if current_time_step % graph_time == 0:
            graph(y_substrate, x_steps, vegf, pedf, fib, pro, x_vector, y_vector, workspace, current_time_step,
                  total_number_time_steps, total_time)

    graph(y_substrate, x_steps, vegf, pedf, fib, pro, x_vector, y_vector, workspace, current_time_step,
          total_number_time_steps, total_time)

    return
