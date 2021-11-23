# Description
# proliferation takes the divide and death probabilities, generates a random number, and then decides if the EC 
# divides, dies, or dies nothing. Then it carries out the necessary action.


# Imports
from random import random


# Function
def proliferation(death_prob, death_time, cell, current_time_step, occupied, y, x, workspace, cell_line, file_deaths,
                divide_prob, number_of_cells, max_cells_allowed, x_position, y_position, total_number_time_steps,
                birth_time, divide_time, file_divisions):


    random_prob = random()  # Random number for dividing or dying

    # Find out if it died and kill it
    if random_prob < death_prob:
        death_time[cell] = current_time_step
        occupied[y][x] -= 1
        workspace[y][x] = cell_line[cell]
        file_deaths.write("\n\nTime: " + str(current_time_step) + "\n")
        file_deaths.write("Cell: " + str(cell) + "\n")
        file_deaths.write("Probability of Death: " + str(death_prob) + "\n")
        file_deaths.write("Random Number: " + str(random_prob))

    # Find out if cell divided and split it
    elif random_prob < death_prob + divide_prob:

        # Don't allow any more than the max amount of cells
        if number_of_cells >= max_cells_allowed:
            print("The max cell limit has been reached, no more divisions")

        # Divide the cell and create a new one
        else:
            x_position[number_of_cells][current_time_step + 1] = x
            y_position[number_of_cells][current_time_step + 1] = y
            occupied[y][x] += 1
            cell_line[number_of_cells] = cell_line[cell]
            death_time[number_of_cells] = total_number_time_steps - 1
            birth_time[number_of_cells] = current_time_step
            divide_time[cell] = current_time_step
            divide_time[number_of_cells] = current_time_step
            number_of_cells += 1
            file_divisions.write("\n\nTime: " + str(current_time_step) + "\n")
            file_divisions.write("Cell: " + str(cell) + "\n")
            file_divisions.write("Probability of Division: " + str(divide_prob) + "\n")
            file_divisions.write("Random Number: " + str(random_prob))

    return
