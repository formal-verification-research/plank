# Description
# proliferation takes the divide and death probabilities, generates a random number, and then decides if the EC 
# divides, dies, or dies nothing. Then it carries out the necessary action.


# Imports
from random import random
from math import exp
from the_vault import K6, K18, K20, K25, K26, M1
from heaviside import heaviside


# Function
def proliferation(death_time, cell, current_time_step, occupied, y, x, workspace, cell_line, file_events,
                  number_of_cells, max_cells_allowed, x_position, y_position, total_number_time_steps, birth_time,
                  divide_time, pro, pedf, fib, pro_old, pedf_old, fib_old, x_steps, y_steps, k, child):
    
    pro0 = 0  # pro values at current time step
    pro1 = 0  # pro values at current time step - 1
    surrounding_points = 0

    if x > 0:  # Left
        pro0 += pro[2 * y][x - 1] / (1 + K6 * fib[2 * y][x - 1])
        pro1 += pro_old[2 * y][x - 1] / (1 + K6 * fib_old[2 * y][x - 1])
        surrounding_points += 1
    if x < x_steps - 1:  # Right
        pro0 += pro[2 * y][x] / (1 + K6 * fib[2 * y][x])
        pro1 += pro_old[2 * y][x] / (1 + K6 * fib_old[2 * y][x])
        surrounding_points += 1
    if y > 0:  # Up
        pro0 += pro[2 * y - 1][x] / (1 + K6 * fib[2 * y - 1][x])
        pro1 += pro_old[2 * y - 1][x] / (1 + K6 * fib_old[2 * y - 1][x])
        surrounding_points += 1
    if y < y_steps - 1:  # Down
        pro0 += pro[2 * y + 1][x - 1] / (1 + K6 * fib[2 * y + 1][x - 1])
        pro1 += pro_old[2 * y + 1][x - 1] / (1 + K6 * fib_old[2 * y + 1][x - 1])
        surrounding_points += 1

    pro0 = pro0 / surrounding_points
    pro1 = pro1 / surrounding_points
    G = K25 * (exp(-K26 * pro0 ** M1) * (1 - K26 * M1 * pro0 ** M1)) / (1 + K25 * pro0 * exp(-K26 * pro0 ** M1))  # G
    pro_dependent = G * (pro0 - pro1) / k  # pro dependant term for divide and death probabilities

    if pro_dependent >= 0:  # Divide and death probabilities
        divide_prob = (k * K18 + G * (pro0 - pro1)) * heaviside(current_time_step - divide_time[cell] - child)
        death_prob = k * K20
    else:
        divide_prob = k * K18 * heaviside(current_time_step - divide_time[cell] - child)
        death_prob = k * K20 - G * (pro0 - pro1)

    random_prob = random()  # Random number for dividing or dying

    if random_prob < death_prob:  # Find out if it died and kill it
        death_time[cell] = current_time_step
        occupied[y][x] -= 1
        workspace[y][x] = cell_line[cell]
        file_events.write("\n\nTime: " + str(current_time_step) + "\n")
        file_events.write("Cell: " + str(cell) + "\n")
        file_events.write("Probability of Death: " + str(death_prob) + "\n")
        file_events.write("Random Number: " + str(random_prob))

    elif random_prob < death_prob + divide_prob:  # Find out if cell divided and split it

        if number_of_cells >= max_cells_allowed:  # Don't allow any more than the max amount of cells
            print("The max cell limit has been reached, no more divisions")

        else:  # Divide the cell and create a new one
            x_position[number_of_cells][current_time_step + 1] = x
            y_position[number_of_cells][current_time_step + 1] = y
            occupied[y][x] += 1
            cell_line[number_of_cells] = cell_line[cell]
            death_time[number_of_cells] = total_number_time_steps - 1
            birth_time[number_of_cells] = current_time_step
            divide_time[cell] = current_time_step
            divide_time[number_of_cells] = current_time_step
            number_of_cells += 1
            file_events.write("\n\nTime: " + str(current_time_step) + "\n")
            file_events.write("Cell: " + str(cell) + "\n")
            file_events.write("Probability of Division: " + str(divide_prob) + "\n")
            file_events.write("Random Number: " + str(random_prob))

    return
