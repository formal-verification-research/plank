# Description
# update_vegf updates the vegf substrate for each new time step


# Imports
from numpy import zeros
from math import cos
from math import pi
from x_coordinate import x_coordinate
from the_vault import K1, K2, K35, K21, M0, K33, RELAX1


# Function
def update_vegf(y_substrate, x_steps, density_scale, occupied_old, vegf, vegf_old, k, tolerance, h, x_length):

    # Initialize v and v_old: substrate matrices used to iterate
    v = zeros((y_substrate, x_steps))
    v_old = zeros((y_substrate, x_steps))
    for y in range(1, y_substrate, 1):
        if y % 2 == 0:
            for x in range(x_steps - 1):
                v[y][x] = vegf[y][x]
        else:
            for x in range(x_steps):
                v[y][x] = vegf[y][x]

    # Capillary
    for x in range(x_steps - 1):
        density = density_scale * (occupied_old[0][x] + occupied_old[0][x + 1]) / 2  # Ave density right and left
        wall_vegf = (vegf[1][x] + vegf[1][x + 1]) / 2  # Average VEGF concentration on capillary wall at time j
        vegf_difference = wall_vegf - vegf[0][x]
        if vegf_difference < 0:
            vegf_difference = 0
        vegf_old[0][x] = vegf[0][x]
        vegf[0][x] = vegf[0][x] + k * (-K1 * vegf[0][x] * density
                                       / (1 + vegf[0][x]) + K2 * vegf_difference)  # plank eq 46
        if vegf[0][x] < 0:
            vegf[0][x] = 0

    # ECM
    in_tol = 0
    while in_tol == 0:
        in_tol = 1

        # Iterate v into v_old
        for y in range(1, y_substrate, 1):
            if y % 2 == 0:
                for x in range(x_steps - 1):
                    v_old[y][x] = v[y][x]
            else:
                for x in range(x_steps):
                    v_old[y][x] = v[y][x]

        # Y = Max, Source
        v[y_substrate - 1][0] = v[y_substrate - 1][1]  # plank pg 179 eq 70
        for x in range(1, x_steps - 2, 1):
            v[y_substrate - 1][x] = K35 * h \
                                    * ((1 - cos(2 * pi * x_coordinate(x, y_substrate - 1, x_steps, x_length))) ** M0) \
                                    + v[y_substrate - 3][x]  # plank pg 179 eq 65
        v[y_substrate - 1][x_steps - 2] = v[y_substrate - 1][x_steps - 3]  # plank pg 179 eq 70

        # Y = Max - 1, Source
        v[y_substrate - 2][0] = v[y_substrate - 2][1]  # plank pg 179 eq 70
        for x in range(1, x_steps - 1, 1):
            v[y_substrate - 2][x] = K35 * h \
                                    * ((1 - cos(2 * pi * x_coordinate(x, y_substrate - 2, x_steps, x_length))) ** M0) \
                                    + v[y_substrate - 4][x]  # plank pg 179 eq 65
        v[y_substrate - 2][x_steps - 1] = v[y_substrate - 2][x_steps - 2]  # plank pg 179 eq 70

        # Interior
        for y in range(y_substrate - 3, 2, -1):
            v[y][0] = v[y][1]  # plank pg 179 eq 70

            if y % 2 == 0:  # If row is even number of substrate mesh points in x is x_steps-1
                for x in range(1, x_steps - 2, 1):  # minus 2 because last mesh point has a boundary condition
                    density = density_scale * ((y_substrate / 2) - 1) \
                              * (occupied_old[y//2][x] + occupied_old[y // 2][x + 1]) / 2  # Av density right and left
                    v[y][x] = RELAX1 / (h * h + 2 * K21 * k) \
                              * (0.5 * K21 * k * (v[y][x + 1] + v[y][x - 1] + v[y + 2][x] + v[y - 2][x]
                                                  + vegf[y][x + 1] + vegf[y][x - 1] + vegf[y + 2][x] + vegf[y - 2][x])
                                 + (h * h - 2 * K21 * k - h * h * k * K1 * density / (1 + vegf[y][x])) * vegf[y][x]) \
                              + (1 - RELAX1) * v[y][x]  # plank pg 178 eq 53 crank-nicolson approximation
                v[y][x_steps - 2] = v[y][x_steps - 3]  # plank pg 179 eq 70

            else:
                for x in range(1, x_steps - 1, 1):  # If row is odd number of substrate mesh points in x is nn
                    density = density_scale * ((y_substrate / 2) - 1) \
                              * (occupied_old[(y - 1) // 2][x] + occupied_old[(y + 1) // 2][x]) / 2  # Ave den above and below
                    v[y][x] = RELAX1 / (h * h + 2 * K21 * k) \
                              * (0.5 * K21 * k * (v[y][x + 1] + v[y][x - 1] + v[y + 2][x] + v[y - 2][x]
                                                  + vegf[y][x + 1] + vegf[y][x - 1] + vegf[y + 2][x] + vegf[y - 2][x])
                                 + (h * h - 2 * K21 * k - h * h * k * K1 * density / (1 + vegf[y][x])) * vegf[y][x]) \
                              + (1 - RELAX1) * v[y][x]
                v[y][x_steps - 1] = v[y][x_steps - 2]  # plank pg 179 eq 70

        # Y = 2: capillary wall
        v[2][0] = v[2][1]  # plank pg 179 eq 70
        for x in range(1, x_steps - 2, 1):
            v[2][x] = 1 / (K33 + h) * (K33 * v[4][x] + h * vegf_old[0][x])  # plank pg 179, 137 eq 61
        v[2][x_steps - 2] = v[2][x_steps - 3]  # plank pg 179 eq 70

        # Y = 1: capillary wall
        v[1][0] = v[1][1]
        for x in range(1, x_steps - 1, 1):
            v[1][x] = 1 / (K33 + h) * (K33 * v[3][x] + h * (vegf_old[0][x - 1] + vegf_old[0][x]) / 2)  # Ave
        v[1][x_steps - 1] = v[1][x_steps - 2]  # plank pg 179 eq 70

        # Check to make sure each point meets the tolerance
        for y in range(1, y_substrate, 1):
            if y % 2 == 0:
                for x in range(x_steps - 1):
                    if v[y][x] - v_old[y][x] > tolerance or v[y][x] - v_old[y][x] < -tolerance:
                        in_tol = 0
            else:
                for x in range(x_steps):
                    if v[y][x] - v_old[y][x] > tolerance or v[y][x] - v_old[y][x] < -tolerance:
                        in_tol = 0

    # Cycle through VEGF mesh points and set VEGF at time step j+1
    for y in range(1, y_substrate, 1):  # If y is even: number of substrate mesh points in x is x_steps - 1
        if y % 2 == 0:
            for x in range(x_steps - 1):
                vegf_old[y][x] = vegf[y][x]
                vegf[y][x] = v[y][x]
                if vegf[y][x] < 0:
                    vegf[y][x] = 0
        else:  # If y is odd: number of substrate mesh points in x is x_steps
            for x in range(x_steps):
                vegf_old[y][x] = vegf[y][x]
                vegf[y][x] = v[y][x]
                if vegf[y][x] < 0:
                    vegf[y][x] = 0

    return vegf, vegf_old
