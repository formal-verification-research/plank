'''
Description:
# The startup file takes the parameters from the parameter vault and builds the model used for the simulation.
'''


# Imports
import time
from numpy import zeros
from numpy import ones
from numpy import linspace
from parameter_vault import x_steps, simulation_time, DP, L, max_cells_allowed, number_of_cells, threshold_perc, \
    division, time_step_duration


# Function
def startup():

    # Define Important Non-dimensionalized Terms
    start_time = time.time()  # Time that the code started running
    x_length = 1  # L
    y_length = 0.5  # Depth, 1/2 of L
    y_steps = int(x_steps * y_length + 0.5)  # y_steps is used for the ec arrays and has 1/2 the steps = 101
    y_substrate = x_steps  # y_substrate is used for substrate arrays, same amount of nodes as x_steps = 201
    file_events = open("EC_Events.txt", "w")  # Stores info on EC divisions, deaths, and leaving the parent blood vessel
    total_time = simulation_time * DP / (L ** 2)  # Time
    total_number_time_steps = int(simulation_time * 3600 / time_step_duration)  # The time steps used in the simulation
    k = total_time / total_number_time_steps  # The time step interval
    h = x_length / x_steps  # The distance between mesh points
    lam = (x_length ** 2) / (h ** 2)  # Size ratio between simulation length and step size (Plank paper pg 150).
    density_cap = x_steps / number_of_cells  # Density of cells in the capillary
    density_ecm = x_steps * (y_steps - 1) / number_of_cells  # Density of cells in the ECM
    threshold = threshold_perc / 100  # Threshold percentage as an amount of Fibronectin
    child = division * 3600 / time_step_duration  # Minimum amount of time permitted between EC divisions in time steps

    # Create vectors used for 3D graphs
    x_vector = []
    x_even_vector = []
    x_odd_vector = []
    y_vector = []
    for x in linspace(0.5, x_steps - 1.5, num=x_steps - 1):
        x_even_vector.append(x)
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

    # Create EC and Substrate Arrays
    x_position = zeros((max_cells_allowed, total_number_time_steps))  # X position vector
    y_position = zeros((max_cells_allowed, total_number_time_steps))  # Y position vector
    death_time = zeros(max_cells_allowed)  # Death array
    birth_time = zeros(max_cells_allowed)  # Birth array
    divide_time = zeros(max_cells_allowed)  # Divide array
    vegf = zeros((y_substrate, x_steps))  # VEGF array
    pedf = zeros((y_substrate, x_steps))  # PEDF array
    pro = zeros((y_substrate, x_steps))  # Pro array
    fib = ones((y_substrate, x_steps))  # Fib array
    vegf_old = zeros((y_substrate, x_steps))  # Old VEGF array
    pedf_old = zeros((y_substrate, x_steps))  # Old PEDF array
    pro_old = zeros((y_substrate, x_steps))  # Old Pro array
    fib_old = ones((y_substrate, x_steps))  # Old Fib array
    model = zeros((y_steps, x_steps))  # EC array for graphing
    ec = zeros((y_steps, x_steps))  # EC array
    ec_old = zeros((y_steps, x_steps))  # Old EC array

    # Trying to create different colors for the different starting EC
    cell_lineage = zeros(max_cells_allowed)
    cl = 10
    for i in range(number_of_cells):
        cell_lineage[i] = cl
        cl += 10

    # Track the EC for use in data outputs such as speed and distance
    cell_tracker = [[[], [], []] for i in range(max_cells_allowed)]

    # Seed initial EC, spaced out evenly along the length of the parent blood vessel
    cell_number = 0  # Used to place all the EC and set up the matrices
    for x in range(int(0.5 * x_steps / number_of_cells), x_steps, int(x_steps / number_of_cells)):
        x_position[cell_number][0] = x
        y_position[cell_number][0] = 0
        ec[0][x] += 1
        death_time[cell_number] = total_number_time_steps - 1  # Set the death time of all EC to the max time
        cell_number += 1

    return start_time, x_length, y_length, y_steps, y_substrate, file_events, total_time, total_number_time_steps, k, \
           h, lam, x_vector, y_vector, x_position, y_position, death_time, birth_time, divide_time, vegf, pedf, pro, \
           fib, vegf_old, pedf_old, pro_old, fib_old, model, ec, ec_old, density_cap, density_ecm, cell_lineage, \
           cell_tracker, cell_number, threshold, child
