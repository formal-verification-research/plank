"""
Description
The update_pedf file updates the PEDF substrate array for the next time step
"""


# Imports
from numpy import zeros
from math import cos
from math import pi
from x_coordinate import x_coordinate
from parameter_vault import K36, K37, K38, K39, M0, K40, RELAX1, x_steps, tolerance


# Function
def update_pedf(y_substrate, density_cap, density_ecm, ec_old, pedf, pedf_old, k, h, x_length, current_time_step, total_number_time_steps):

    # Create new PEDF arrays to help iterate in the ECM iterative equations
    p = zeros((y_substrate, x_steps))
    p_old = zeros((y_substrate, x_steps))
    for y in range(1, y_substrate):
        if y % 2 == 0:
            for x in range(x_steps-1):
                p[y][x] = pedf[y][x]
        else:
            for x in range(x_steps):
                p[y][x] = pedf[y][x]

    # Update the PEDF inside the parent blood vessel
    for x in range(x_steps-1):
        pedf_old[0][x] = pedf[0][x]
        density = density_cap * (ec_old[0][x] + ec_old[0][x+1]) / 2
        wall_pedf = (pedf[1][x] + pedf[1][x+1]) / 2
        pedf_difference = wall_pedf - pedf[0][x]
        if pedf_difference < 0:
            pedf_difference = 0
        pedf[0][x] = pedf[0][x] + k * (-K36 * pedf[0][x] * density / (1 + pedf[0][x]) + K37 * pedf_difference)
        if pedf[0][x] < 0:
            pedf[0][x] = 0

    # Update the PEDF in the blood vessel wall and in the ECM
    tol = 0
    while tol == 0:
        tol = 1

        # Take the current iteration and place it into the PEDF array
        for y in range(1, y_substrate):
            if y % 2 == 0:
                for x in range(x_steps-1):
                    p_old[y][x] = p[y][x]
            else:
                for x in range(x_steps):
                    p_old[y][x] = p[y][x]

        # Update the PEDF source, which is the RPE layer and just outside
        for x in range(1, x_steps-2):
            p[y_substrate-1][x] = K38 * h \
                                    * ((1 - cos(2 * pi * x_coordinate(x, y_substrate - 1, x_steps, x_length))) ** M0) \
                                    + p[y_substrate-3][x]
        p[y_substrate-1][0] = p[y_substrate-1][1]
        p[y_substrate-1][x_steps-2] = p[y_substrate-1][x_steps-3]
        for x in range(1, x_steps-1):
            p[y_substrate-2][x] = K38 * h \
                                    * ((1 - cos(2 * pi * x_coordinate(x, y_substrate - 2, x_steps, x_length))) ** M0) \
                                    + p[y_substrate-4][x]
        p[y_substrate-2][0] = p[y_substrate-2][1]
        p[y_substrate-2][x_steps-1] = p[y_substrate-2][x_steps-2]

        # Update the PEDF interior points (Crank-Nicolson Approximation)
        for y in range(y_substrate-3, 2, -1):
            if y % 2 == 0:
                for x in range(1, x_steps-2):
                    density = density_ecm * (ec_old[y//2][x] + ec_old[y//2][x+1]) / 2
                    p[y][x] = RELAX1 / (h * h + 2 * K39 * k) \
                              * (0.5 * K39 * k * (p[y][x+1] + p[y][x-1] + p[y+2][x] + p[y-2][x]
                                    + pedf[y][x+1] + pedf[y][x-1] + pedf[y+2][x] + pedf[y-2][x])
                                 + (h * h - 2 * K39 * k - h * h * k * K36 * density / (1 + pedf[y][x]))
                                 * pedf[y][x]) + (1 - RELAX1) * p[y][x]
                p[y][0] = p[y][1]
                p[y][x_steps-2] = p[y][x_steps-3]
            else:
                for x in range(1, x_steps-1):
                    density = density_ecm * (ec_old[(y-1)//2][x] + ec_old[(y+1)//2][x]) / 2
                    p[y][x] = RELAX1 / (h * h + 2 * K39 * k) \
                              * (0.5 * K39 * k * (p[y][x+1] + p[y][x-1] + p[y+2][x] + p[y-2][x]
                                    + pedf[y][x+1] + pedf[y][x-1] + pedf[y+2][x] + pedf[y-2][x])
                                 + (h * h - 2 * K39 * k - h * h * k * K36 * density / (1 + pedf[y][x]))
                                 * pedf[y][x]) \
                              + (1 - RELAX1) * p[y][x]
                p[y][0] = p[y][1]
                p[y][x_steps-1] = p[y][x_steps-2]

        # Update the PEDF on the parent blood vessel wall
        for x in range(1, x_steps-2):
            p[2][x] = 1 / (K40 + h) * (K40 * p[4][x] + h * pedf_old[0][x])
        p[2][0] = p[2][1]
        p[2][x_steps-2] = p[2][x_steps-3]
        for x in range(1, x_steps-1):
            p[1][x] = 1 / (K40 + h) * (K40 * p[3][x] + h * (pedf_old[0][x] + pedf_old[0][x-1]) / 2)
        p[1][0] = p[1][1]
        p[1][x_steps-1] = p[1][x_steps-2]

        # Check to make sure each point meets the tolerance before the while loop breaks
        for y in range(1, y_substrate):
            if y % 2 == 0:
                for x in range(x_steps-1):
                    if p[y][x] - p_old[y][x] > tolerance or p[y][x] - p_old[y][x] < -tolerance:
                        tol = 0
            else:
                for x in range(x_steps):
                    if p[y][x] - p_old[y][x] > tolerance or p[y][x] - p_old[y][x] < -tolerance:
                        tol = 0

    # Cycle through PEDF mesh points and set PEDF for the future time step and previous time step
    for y in range(1, y_substrate):
        if y % 2 == 0:
            for x in range(x_steps-1):
                pedf_old[y][x] = pedf[y][x]
                pedf[y][x] = p[y][x]
                if current_time_step / total_number_time_steps > .05:
                    pedf_old[y][x] = pedf[y][x] * .99
                    pedf[y][x] = p[y][x] * .99

                if pedf[y][x] < 0:
                    pedf[y][x] = 0
        else:
            for x in range(x_steps):
                pedf_old[y][x] = pedf[y][x]
                pedf[y][x] = p[y][x]
                if current_time_step / total_number_time_steps > .05:
                    pedf_old[y][x] = pedf[y][x] * .99
                    pedf[y][x] = p[y][x] * .99
                if pedf[y][x] < 0:
                    pedf[y][x] = 0

    return pedf, pedf_old
