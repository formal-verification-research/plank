# Description
# start takes the information from the main file and creates the parameters, vectors, and arrays needed for the
# simulation. It creates the EC and substrate matrices and places the EC inside the parent capillary.


# Imports
import time
from numpy import zeros
from numpy import ones
from numpy import linspace


# Function
def terms(x_steps, total_time, total_number_time_steps):
    start_time = time.time()  # Lock in the start time of the program
    x_length = 1  # Nondimensionalize the length
    y_length = 0.5  # Height is 1/2 the length
    y_steps = int(x_steps * (y_length / x_length) + 0.5)  # Half the steps, used for the EC
    y_substrate = y_steps * 2 - 1  # Full steps, used for substrates
    file_events = open("EC_Events.txt", "w")  # To write deaths and divisions to
    k = total_time / total_number_time_steps  # The time step interval
    h = x_length / x_steps  # Create h, the distance between mesh points variable
    lam = (x_length ** 2) / (h ** 2)  # Ratio between length and step size. Plank paper pg 150
    x_vector = []  # used to create 3D graphs
    y_vector = []  # used to create 3D graphs
    x_even_vector = []  # Place even rows

    for x in linspace(0.5, x_steps - 1.5, num=x_steps - 1):
        x_even_vector.append(x)
    x_odd_vector = []  # Place odd rows
    for x in linspace(0, x_steps - 1, num=x_steps):
        x_odd_vector.append(x)
    for y in range(y_substrate):  # Place into vectors
        if y % 2 == 0:
            for i in x_even_vector:
                x_vector.append(i)
                y_vector.append(y)
        else:
            for j in x_odd_vector:
                x_vector.append(j)
                y_vector.append(y)

    return start_time, x_length, y_length, y_steps, y_substrate, file_events, k, h, lam, x_vector, y_vector


# Function
def arrays(max_cells_allowed, total_number_time_steps, y_substrate, x_steps, y_steps):
    x_position = zeros((max_cells_allowed, total_number_time_steps), dtype=int)  # x position vector
    y_position = zeros((max_cells_allowed, total_number_time_steps), dtype=int)  # y position vector
    death_time = zeros(max_cells_allowed, dtype=int)  # Death array
    birth_time = zeros(max_cells_allowed, dtype=int)  # Birth array
    divide_time = zeros(max_cells_allowed, dtype=int)  # Divide array
    vegf = zeros((y_substrate, x_steps))  # vegf array
    pedf = zeros((y_substrate, x_steps))  # pedf array
    pro = zeros((y_substrate, x_steps))  # pro array
    fib = ones((y_substrate, x_steps))  # fib array
    vegf_old = zeros((y_substrate, x_steps))  # old vegf array
    pedf_old = zeros((y_substrate, x_steps))  # old pedf array
    pro_old = zeros((y_substrate, x_steps))  # old pro array
    fib_old = ones((y_substrate, x_steps))  # old fib array
    workspace = zeros((y_steps, x_steps))  # Create the EC matrix

    return x_position, y_position, death_time, birth_time, divide_time, vegf, pedf, pro, fib, vegf_old, pedf_old, \
           pro_old, fib_old, workspace


# Function
def parent_vessel(y_steps, x_steps, number_of_cells, max_cells_allowed, x_position, y_position, death_time,
                  total_number_time_steps):
    occupied = zeros((y_steps, x_steps), dtype=int)  # To tell if that space has an EC in it or not
    occupied_old = zeros((y_steps, x_steps), dtype=int)  # To tell if that space has an EC in it or not
    density_scale = x_steps / number_of_cells  # Density of cells in the capillary

    cell_line = zeros(max_cells_allowed)  # Try and create different cell line colors
    cl = 20
    for i in range(number_of_cells):
        cell_line[i] = cl
        cl += 15

    cell_tracker = [[[], [], []] for i in range(max_cells_allowed)]  # Create the cell tracking vector
    cell_number = 0  # Used to place all the EC and set up the matrices

    # Seed initial EC, space out evenly along the length
    for x in range(int(0.5 * x_steps / number_of_cells), x_steps, int(x_steps / number_of_cells)):
        x_position[cell_number][0] = x
        y_position[cell_number][0] = 0
        occupied[0][x] += 1

        # Set the death time of all EC to the max time
        death_time[cell_number] = total_number_time_steps - 1
        cell_number += 1

    return occupied, occupied_old, density_scale, cell_line, cell_tracker, cell_number
