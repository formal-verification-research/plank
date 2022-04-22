"""
Description
The update_vegf file updates the VEGF substrate array for the next time step
"""


# Imports
from numpy import zeros
from math import cos
from math import pi
from x_coordinate import x_coordinate
from parameter_vault import K1, K2, K35, K21, M0, K33, RELAX1, x_steps, tolerance


# Function
def update_vegf(y_substrate, density_cap, density_ecm, ec_old, vegf, vegf_old, k, h, x_length):

    # Update the VEGF inside the parent blood vessel using Plank Eq 46
    for x in range(x_steps-1):
        vegf_old[0][x] = vegf[0][x]
        density = density_cap * (ec_old[0][x] + ec_old[0][x+1]) / 2
        vegf[0][x] = vegf_old[0][x] + k * (-K1 * vegf_old[0][x] * density / (1 + vegf_old[0][x]) + K2
                                           * ((vegf_old[1][x] + vegf_old[1][x+1]) / 2 - vegf_old[0][x]))
        if vegf[0][x] < 0:
            vegf[0][x] = 0

    # Set VEGF for the previous time step
    for y in range(1, y_substrate):
        if y % 2 == 0:
            for x in range(x_steps-1):
                vegf_old[y][x] = vegf[y][x]
        else:
            for x in range(x_steps):
                vegf_old[y][x] = vegf[y][x]

    # Create a VEGF array to help iterate in the ECM equations
    v = zeros((y_substrate, x_steps))
    for y in range(1, y_substrate):
        if y % 2 == 0:
            for x in range(x_steps-1):
                v[y][x] = vegf[y][x]
        else:
            for x in range(x_steps):
                v[y][x] = vegf[y][x]

    # Update the VEGF in the blood vessel wall and in the ECM using iteration loop
    tol = 0
    while tol == 0:
        tol = 1

        # Take the current iteration and place it into the VEGF array
        for y in range(1, y_substrate):
            if y % 2 == 0:
                for x in range(x_steps-1):
                    vegf[y][x] = v[y][x]
            else:
                for x in range(x_steps):
                    vegf[y][x] = v[y][x]

        # Update the VEGF source, which is the RPE layer and just outside, using Plank Eq 65 and 70 on Pg 179
        for x in range(x_steps-1):
            v[y_substrate-1][x] = K35 * h \
                                    * ((1 - cos(2 * pi * x_coordinate(x, y_substrate - 1, x_steps, x_length))) ** M0) \
                                    + v[y_substrate-3][x]
        for x in range(x_steps):
            v[y_substrate-2][x] = K35 * h \
                                    * ((1 - cos(2 * pi * x_coordinate(x, y_substrate - 2, x_steps, x_length))) ** M0) \
                                    + v[y_substrate-4][x]

        # Update the VEGF interior points using Plank Eq 53, 70 Pgs 178-179 (Crank-Nicolson Approximation)
        for y in range(y_substrate-3, 2, -1):
            if y % 2 == 0:
                for x in range(1, x_steps-2):
                    density = density_ecm * (ec_old[y//2][x] + ec_old[y//2][x+1]) / 2
                    v[y][x] = RELAX1 \
                              * k / (h * h + 2 * K21 * k) \
                              * ((K21 / 2)
                                 * (v[y][x+1] + v[y][x-1] + v[y+2][x] + v[y-2][x]
                                    + vegf_old[y][x+1] + vegf_old[y][x-1] + vegf_old[y+2][x] + vegf_old[y-2][x])
                                 + (h * h / k - 2 * K21 - h * h * K1 * density / (1 + vegf_old[y][x]))
                                 * vegf_old[y][x]) \
                              + (1 - RELAX1) * v[y][x]
                v[y][0] = v[y][1]
                v[y][x_steps-2] = v[y][x_steps-3]
            else:
                for x in range(1, x_steps-1):
                    density = density_ecm * (ec_old[(y-1)//2][x] + ec_old[(y+1)//2][x]) / 2
                    v[y][x] = RELAX1 \
                              * k / (h * h + 2 * K21 * k) \
                              * ((K21 / 2)
                                 * (v[y][x+1] + v[y][x-1] + v[y+2][x] + v[y-2][x]
                                    + vegf_old[y][x+1] + vegf_old[y][x-1] + vegf_old[y+2][x] + vegf_old[y-2][x])
                                 + (h * h / k - 2 * K21 - h * h * K1 * density / (1 + vegf_old[y][x]))
                                 * vegf_old[y][x]) \
                              + (1 - RELAX1) * v[y][x]
                v[y][0] = v[y][1]
                v[y][x_steps-1] = v[y][x_steps-2]

        # Update the VEGF on the parent blood vessel wall using Plank Eq 61, 70 Pgs 137, 179
        for x in range(1, x_steps-2):
            v[2][x] = 1 / (K33 + h) * (K33 * v[4][x] + h * v[0][x])
            v[2][0] = v[2][1]
            v[2][x_steps-2] = v[2][x_steps-3]
        for x in range(1, x_steps-1):
            v[1][x] = 1 / (K33 + h) * (K33 * v[3][x] + h * (v[0][x-1] + v[0][x]) / 2)
            v[1][0] = v[1][1]
            v[1][x_steps-1] = v[1][x_steps-2]

        # Check to make sure each point meets the tolerance before the while loop breaks
        for y in range(1, y_substrate):
            if y % 2 == 0:
                for x in range(x_steps-1):
                    if v[y][x] - vegf[y][x] > tolerance or v[y][x] - vegf[y][x] < -tolerance:
                        tol = 0
            else:
                for x in range(x_steps):
                    if v[y][x] - vegf[y][x] > tolerance or v[y][x] - vegf[y][x] < -tolerance:
                        tol = 0

    # Cycle through VEGF mesh points and set VEGF for the future time step
    for y in range(1, y_substrate):
        if y % 2 == 0:
            for x in range(x_steps-1):
                vegf[y][x] = v[y][x]
                if vegf[y][x] < 0:
                    vegf[y][x] = 0
        else:
            for x in range(x_steps):
                vegf[y][x] = v[y][x]
                if vegf[y][x] < 0:
                    vegf[y][x] = 0

    return vegf, vegf_old
