"""
Description
The simulation file takes the model created by the startup file and runs the simulation until completion.
"""


# Imports
from tip_cell import tip_cell


# Function
def simulation(start_time, x_length, y_length, y_steps, y_substrate, file_events, total_time, total_number_time_steps, k, \
    h, lam, x_vector, y_vector, x_position, y_position, death_time, birth_time, divide_time, vegf, pedf, pro, \
    fib, vegf_old, pedf_old, pro_old, fib_old, model, ec, ec_old, density_cap, density_ecm, cell_lineage, \
    cell_tracker, cell_number, number_of_cells):

    # Start the 'for-loop' that will take the simulation through the time steps
    for current_time_step in range(total_number_time_steps - 1):

        # Start the 'for-loop' that will cycle through each EC for each time step
        for cell in range(number_of_cells):

            # Use the tip_cell function to run the processes that happen to a single EC for a single time step
            tip_cell()







    # # Cycle through time steps
    # for current_time_step in range(total_number_time_steps - 1):

        # # Copy occupied into occupied_old
        # for x in range(x_steps):
        #     for y in range(y_steps):
        #         occupied_old[y][x] = occupied[y][x]

        # # Cycle through cells
        # for cell in range(number_of_cells):

            # # At first assume no movement
            # x = x_position[cell][current_time_step]
            # y = y_position[cell][current_time_step]
            # x_position[cell][current_time_step + 1] = x
            # y_position[cell][current_time_step + 1] = y

            # # If cell has left the capillary
            # if death_time[cell] == total_number_time_steps - 1 and y > 0:

                # # Add the time, x, y to the cell tracking vector
                # cell_tracker[cell][0].append(current_time_step)
                # cell_tracker[cell][1].append(x)
                # cell_tracker[cell][2].append(y)

                # # The cell dies/leaves simulation if it reaches the tumour
                # if y == y_steps - 1:
                #     death_time[cell] = current_time_step
                #     occupied[y][x] -= 1

                # number_of_cells = \
                #     proliferation(death_time, cell, current_time_step, occupied, y, x, workspace, cell_line,
                #                   file_events, number_of_cells, max_cells_allowed, x_position, y_position,
                #                   total_number_time_steps, birth_time, divide_time, pro, pedf, fib, pro_old, pedf_old,
                #                   fib_old, x_steps, y_steps, k, child)

            # # Determine if the cell moves and where
            # if death_time[cell] == total_number_time_steps - 1:
            #     stay = prob_stay(y, lam, k)
            #     left, T = prob_move(x, y, 0, pro, fib, vegf, x_steps, y_steps, lam, k)
            #     right, T = prob_move(x, y, 1, pro, fib, vegf, x_steps, y_steps, lam, k)
            #     up, T = prob_move(x, y, 2, pro, fib, vegf, x_steps, y_steps, lam, k)
            #     random_num = random()

                # # Check if cell can escape the capillary
                # if y == 0:
                #     if x == 0:
                #         fib_cap = fib[0][0]
                #     elif x == x_steps - 1:
                #         fib_cap = fib[0][x_steps - 2]
                #     else:
                #         fib_cap = (fib[0][x - 1] + fib[0][x]) / 2
                #     if fib_cap < threshold:
                #         random_num = 2
                #         file_events.write("\n\nTime: " + str(current_time_step) + "\n")
                #         file_events.write("Cell: " + str(cell) + "\n")
                #         file_events.write("Left the capillary" + "\n")

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

        vegf, vegf_old = \
            update_vegf(y_substrate, x_steps, density_scale, occupied_old, vegf, vegf_old, k, tolerance, h, x_length)
        pedf, pedf_old = \
            update_pedf(y_substrate, x_steps, density_scale, occupied_old, pedf, pedf_old, k, tolerance, h, x_length)
        fib, fib_old = \
            update_fib(y_substrate, x_steps, density_scale, occupied_old, fib, fib_old, k, pro, tolerance, h)
        pro, pro_old = \
            update_pro(y_substrate, x_steps, density_scale, occupied_old, pro, pro_old, k, vegf_old, pedf_old)

        print("Current Time Step = " + str(current_time_step))

        if current_time_step % graph_time == 0:
            graph(y_substrate, x_steps, vegf, pedf, fib, pro, x_vector, y_vector, workspace, current_time_step,
                  total_number_time_steps, total_time)

    graph(y_substrate, x_steps, vegf, pedf, fib, pro, x_vector, y_vector, workspace, current_time_step,
          total_number_time_steps, total_time)

    return
