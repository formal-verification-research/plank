# Description
# proliferation takes the divide and death probabilities, generates a random number, and then decides if the EC 
# divides, dies, or dies nothing. Then it carries out the necessary action.


# Function
def proliferation(x, y, x_position, y_position, occupied, death_time, birth_time, divide_time, number_of_cells, cell, current_time_step,
             totalNumberOfTimeSteps, divideProbability, fileDivisions, cellLine):

    # Divide the cell and create a new one
    x_position[number_of_cells][current_time_step + 1] = x
    y_position[number_of_cells][current_time_step + 1] = y
    occupied[y][x] += 1
    cellLine[number_of_cells] = cellLine[cell]
    death_time[number_of_cells] = totalNumberOfTimeSteps - 1
    birth_time[number_of_cells] = current_time_step
    divide_time[cell] = current_time_step
    divide_time[number_of_cells] = current_time_step
    number_of_cells += 1

    # Get some data into the divison file
    fileDivisions.write("\n\nTime: " + str(current_time_step) + "\n")
    fileDivisions.write("Cell: " + str(cell) + "\n")
    fileDivisions.write("Probability of Division: " + str(divideProbability) + "\n")

    return
