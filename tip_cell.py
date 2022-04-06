"""
Description:
The tip_cell file runs through all the processes that happen to a single EC during a single time step
"""


# Imports
from proliferation import proliferation
from random import random
from prob_stay import prob_stay
from prob_move import prob_move


# Function
def tip_cell(x_position, cell, current_time_step, y_position, ec_old, ec, death_time,
             total_number_time_steps, cell_tracker, y_steps, pro, pro_old, fib, fib_old, x_steps, k, divide_time,
             child, model, cell_lineage, file_events, number_of_cells, max_cells_allowed, birth_time, lam, vegf,
             threshold):

    # At first assume the EC will not move
    x = x_position[cell][current_time_step]
    y = y_position[cell][current_time_step]
    x_position[cell][current_time_step + 1] = x
    y_position[cell][current_time_step + 1] = y
    ec_old[y][x] = ec[y][x]

    # Perform the following actions if the EC is alive, which means the death time is still set at the end
    if death_time[cell] == total_number_time_steps - 1:

        # Perform the following actions only if the EC is in the ECM, not in the parent blood vessel
        if y > 0:

            # Add the current time step, Y and X EC positions to the cell tracking vector to find speeds and distances
            cell_tracker[cell][0].append(current_time_step)
            cell_tracker[cell][1].append(x)
            cell_tracker[cell][2].append(y)

            # Kill the EC if it reaches the RPE layer so it doesn't start backtracking
            if y == y_steps - 1:
                death_time[cell] = current_time_step
                ec[y][x] -= 1

            # Run the proliferation function to see if the EC will possibly divide or die
            death_time, ec, model, file_events, number_of_cells, x_position, y_position, cell_lineage, death_time, \
            birth_time, divide_time \
                = proliferation(x, pro, pro_old, fib, fib_old, x_steps, y, y_steps, k, current_time_step, divide_time,
                                cell, child, death_time, ec, model, cell_lineage, file_events, number_of_cells,
                                max_cells_allowed, x_position, y_position, total_number_time_steps, birth_time)

        # Find the probability that the EC will stay or move. Generate the random number that decides the EC outcome
        stay = prob_stay(y, lam, k)
        left, T = prob_move(x, y, 0, pro, fib, vegf, x_steps, y_steps, lam, k)
        right, T = prob_move(x, y, 1, pro, fib, vegf, x_steps, y_steps, lam, k)
        up, T = prob_move(x, y, 2, pro, fib, vegf, x_steps, y_steps, lam, k)
        random_num = random()

        # Perform the following action only if the EC is in the parent blood vessel
        if y == 0:

            # Check the Fibronectin levels in the capillary wall to see if the EC can escape the parent blood vessel
            if x == 0:
                fib_cap = fib[0][0]
            elif x == x_steps - 1:
                fib_cap = fib[0][x_steps - 2]
            else:
                fib_cap = (fib[0][x - 1] + fib[0][x]) / 2
            if fib_cap < threshold:
                random_num = 2
                file_events.write("\n\nTime: " + str(current_time_step) + "\n")
                file_events.write("Cell: " + str(cell) + "\n")
                file_events.write("Left the capillary" + "\n")
