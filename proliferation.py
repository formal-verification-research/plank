"""
Description:
The proliferation file samples the average Protease surrounding the EC for the current time step and the previous one.
It then uses these values to find the probability that the EC will divide, die, or neither. These values are places on
a scale from 0-1 and a random number generator decides the outcome, which is then carried out by the function and
returned.
"""


# Imports
from parameter_vault import K6, K18, K20, K25, K26, M1
from math import exp
from heaviside import heaviside
from random import random
from parameter_vault import x_steps, max_cells_allowed


# Function
def proliferation(y_steps, file_events, total_number_time_steps, k, x_position, y_position, death_time, birth_time,
                  divide_time, pro, fib, pro_old, fib_old, model, ec, cell_lineage, number_of_cells, child, cell,
                  current_time_step, x, y, nodes):

    # Find the surrounding pro values at the current and previous time steps
    pro0 = 0  # pro values at current time step
    pro1 = 0  # pro values at previous time step
    surrounding_points = 0  # Nodes surrounding the EC
    # EC location conditions
    if x > 0:  # Left
        pro0 += pro[2*y][x-1] / (1 + K6 * fib[2*y][x-1])
        pro1 += pro_old[2*y][x-1] / (1 + K6 * fib_old[2*y][x-1])
        surrounding_points += 1
    if x < x_steps - 1:  # Right
        pro0 += pro[2*y][x] / (1 + K6 * fib[2*y][x])
        pro1 += pro_old[2*y][x] / (1 + K6 * fib_old[2*y][x])
        surrounding_points += 1
    if y > 0:  # Up
        pro0 += pro[2*y-1][x] / (1 + K6 * fib[2*y-1][x])
        pro1 += pro_old[2*y-1][x] / (1 + K6 * fib_old[2*y-1][x])
        surrounding_points += 1
    if y < y_steps - 1:  # Down
        pro0 += pro[2*y+1][x-1] / (1 + K6 * fib[2*y+1][x-1])
        pro1 += pro_old[2*y+1][x-1] / (1 + K6 * fib_old[2*y+1][x-1])
        surrounding_points += 1
    # Find the average values
    pro0 = pro0 / surrounding_points
    pro1 = pro1 / surrounding_points

    # Find G and GdC/dt, the Protease dependant terms used for calculating probabilities (Plank Paper Pages 153-154)
    G = K25 * (exp(-K26 * pro0 ** M1) * (1 - K26 * M1 * pro0 ** M1)) / (1 + K25 * pro0 * exp(-K26 * pro0 ** M1))
    pro_dependent = G * (pro0 - pro1) / k

    # Find the division and death probabilities based on the Protease values, and generate the random number
    if pro_dependent >= 0:
        divide_prob = (k * K18 + G * (pro0 - pro1)) * heaviside(current_time_step - divide_time[cell] - child)
        death_prob = k * K20
    else:
        divide_prob = k * K18 * heaviside(current_time_step - divide_time[cell] - child)
        death_prob = k * K20 - G * (pro0 - pro1)
    random_prob = random()

    # Carry out the outcome based on the decision made by the random number generator
    # EC death
    if random_prob < death_prob:
        death_time[cell] = current_time_step
        ec[y][x] -= 1
        model[y][x] = cell_lineage[cell]
        file_events.write("\n\nTime: " + str(current_time_step) + "\n")
        file_events.write("Cell: " + str(cell) + " died" + "\n")
        file_events.write("Probability of Death: " + str(death_prob) + "\n")
        file_events.write("Random Number: " + str(random_prob))
    # EC division
    elif random_prob < death_prob + divide_prob:
        # Don't allow any more than the max amount of cells
        if number_of_cells >= max_cells_allowed:
            print("The max cell limit has been reached, no more divisions")
        else:
            x_position[number_of_cells][current_time_step+1] = x
            y_position[number_of_cells][current_time_step+1] = y
            ec[y][x] += 1
            cell_lineage[number_of_cells] = cell_lineage[cell]
            death_time[number_of_cells] = total_number_time_steps - 1
            birth_time[number_of_cells] = current_time_step
            divide_time[cell] = current_time_step
            divide_time[number_of_cells] = current_time_step
            number_of_cells += 1
            file_events.write("\n\nTime: " + str(current_time_step) + "\n")
            file_events.write("Cell: " + str(cell) + " divided" + "\n")
            file_events.write("Probability of Division: " + str(divide_prob) + "\n")
            file_events.write("Random Number: " + str(random_prob))
            nodes += 1

    return number_of_cells, nodes, divide_prob
