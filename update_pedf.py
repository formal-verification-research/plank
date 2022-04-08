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
def update_pedf(y_substrate, density_cap, density_ecm, ec_old, pedf, pedf_old, k, h, x_length):

    # Update the PEDF inside the parent blood vessel
    for x in range(x_steps-1):
        pedf_old[0][x] = pedf[0][x]
        density = density_cap * (ec_old[0][x] + ec_old[0][x+1]) / 2
        pedf[0][x] = pedf_old[0][x] + k * (-K36 * pedf_old[0][x] * density / (1 + pedf_old[0][x]) + K37
                                           * ((pedf_old[1][x] + pedf_old[1][x+1]) / 2 - pedf_old[0][x]))
        if pedf[0][x] < 0:
            pedf[0][x] = 0

    # Create a PEDF array to help iterate in the ECM equations
    p = zeros((y_substrate, x_steps))
    for y in range(1, y_substrate):
        if y % 2 == 0:
            for x in range(x_steps-1):
                p[y][x] = pedf[y][x]
        else:
            for x in range(x_steps):
                p[y][x] = pedf[y][x]

    # Update the PEDF in the blood vessel wall and in the ECM
    tol = 0
    while tol == 0:
        tol = 1

        # Take the current iteration and place it into the PEDF array
        for y in range(1, y_substrate):
            if y % 2 == 0:
                for x in range(x_steps-1):
                    pedf[y][x] = p[y][x]
            else:
                for x in range(x_steps):
                    pedf[y][x] = p[y][x]

        # Update the PEDF source, which is the RPE layer
        for x in range(x_steps-1):
            p[y_substrate-1][x] = K38 * h \
                                    * ((1 - cos(2 * pi * x_coordinate(x, y_substrate - 1, x_steps, x_length))) ** M0) \
                                    + (pedf[y_substrate-2][x] + pedf[y_substrate-2][x+1]) / 2

        # Update the PEDF interior points (Crank-Nicolson Approximation)
        for y in range(y_substrate-2, 1, -1):
            if y % 2 == 0:
                for x in range(1, x_steps-2):
                    density = density_ecm * (ec_old[y//2][x] + ec_old[y//2][x+1]) / 2
                    p[y][x] = RELAX1 \
                              * k / (h * h + 2 * K39 * k) \
                              * ((K39 / 2)
                                 * (pedf[y][x+1] + pedf[y][x-1]
                                    + (pedf[y+1][x] + pedf[y+1][x+1]) / 2 + (pedf[y-1][x] + pedf[y-1][x+1]) / 2
                                    + pedf_old[y][x+1] + pedf_old[y][x-1]
                                    + (pedf_old[y+1][x] + pedf_old[y+1][x+1]) / 2
                                    + (pedf_old[y-1][x] + pedf_old[y-1][x+1]) / 2)
                                 + (h * h / k - 2 * K39 - h * h * K36 * density / (1 + pedf_old[y][x]))
                                 * pedf_old[y][x]) \
                              + (1 - RELAX1) * pedf[y][x]
                p[y][0] = p[y][1]
                p[y][x_steps-2] = p[y][x_steps-3]
            else:
                for x in range(1, x_steps-1):
                    density = density_ecm * (ec_old[(y-1)//2][x] + ec_old[(y+1)//2][x]) / 2
                    p[y][x] = RELAX1 \
                              * k / (h * h + 2 * K39 * k) \
                              * ((K39 / 2)
                                 * (pedf[y][x+1] + pedf[y][x-1]
                                    + (pedf[y+1][x-1] + pedf[y+1][x]) / 2 + (pedf[y-1][x-1] + pedf[y-1][x]) / 2
                                    + pedf_old[y][x+1] + pedf_old[y][x-1]
                                    + (pedf_old[y+1][x-1] + pedf_old[y+1][x]) / 2
                                    + (pedf_old[y-1][x-1] + pedf_old[y-1][x]) / 2)
                                 + (h * h / k - 2 * K39 - h * h * K36 * density / (1 + pedf_old[y][x]))
                                 * pedf_old[y][x]) \
                              + (1 - RELAX1) * pedf[y][x]
                p[y][0] = p[y][1]
                p[y][x_steps-1] = p[y][x_steps-2]

        # Update the PEDF on the parent blood vessel wall
        for x in range(1, x_steps-1):
            p[1][x] = 1 / (K40 + h) * (K40 * (pedf[2][x] + pedf[2][x-1]) / 2 + h * (pedf[0][x] + pedf[0][x-1]) / 2)
            p[1][0] = p[1][1]
            p[1][x_steps-1] = p[1][x_steps-2]

        # Check to make sure each point meets the tolerance before the while loop breaks
        for y in range(1, y_substrate):
            if y % 2 == 0:
                for x in range(x_steps-1):
                    if p[y][x] - pedf[y][x] > tolerance or p[y][x] - pedf[y][x] < -tolerance:
                        tol = 0
            else:
                for x in range(x_steps):
                    if p[y][x] - pedf[y][x] > tolerance or p[y][x] - pedf[y][x] < -tolerance:
                        tol = 0

    # Cycle through PEDF mesh points and set PEDF for the future time step
    for y in range(1, y_substrate):
        if y % 2 == 0:
            for x in range(x_steps-1):
                pedf_old[y][x] = pedf[y][x]
                pedf[y][x] = p[y][x]
                if pedf[y][x] < 0:
                    pedf[y][x] = 0
        else:
            for x in range(x_steps):
                pedf_old[y][x] = pedf[y][x]
                pedf[y][x] = p[y][x]
                if pedf[y][x] < 0:
                    pedf[y][x] = 0

    return pedf, pedf_old
