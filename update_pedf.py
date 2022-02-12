# Description
# update_pedf updates the pedf substrate for each new time step


# Imports
from numpy import zeros
from math import cos
from math import pi
from x_coordinate import x_coordinate
from the_vault import K1, K2, K35, K21, M0, K33, RELAX1


# Function
def update_pedf(y_substrate, x_steps, density_scale, occupied_old, pedf, pedf_old, k, tolerance, h, x_length):

    p = zeros((y_substrate, x_steps))  # Substrate matrix used to iterate

    # Initialize p: value of the previous time step
    for y in range(1, y_substrate, 1):
        if y % 2 == 0:
            for x in range(x_steps - 1):
                p[y][x] = pedf[y][x]
        else:
            for x in range(x_steps):
                p[y][x] = pedf[y][x]

    # Capillary
    for x in range(x_steps - 1):
        density = density_scale * (occupied_old[0][x] + occupied_old[0][x + 1]) / 2  # Ave density right and left
        wall_pedf = (pedf[1][x] + pedf[1][x + 1]) / 2  # Average pedf concentration on capillary wall at time j
        pedf_difference = wall_pedf - pedf[0][x]
        if pedf_difference < 0:
            pedf_difference = 0
        pedf_old[0][x] = pedf[0][x]
        pedf[0][x] = pedf[0][x] + k * (-K1 * pedf[0][x] * density
                                       / (1 + pedf[0][x]) + K2 * pedf_difference)  # plank eq 46
        if pedf[0][x] < 0:
            pedf[0][x] = 0

    # ECM
    in_tol = 0
    while in_tol == 0:
        in_tol = 1

        # Y = Max, includes source

        # X = 0
        p_old = p[y_substrate - 1][0]
        p[y_substrate - 1][0] = p[y_substrate - 1][1]  # plank pg 179 eq 70
        if p[y_substrate - 1][0] - p_old > tolerance or p[y_substrate - 1][0] - p_old < -tolerance:
            in_tol = 0

        # Body
        for x in range(1, x_steps - 2, 1):
            p_old = p[y_substrate - 1][x]
            p[y_substrate - 1][x] = K35 * h \
                                    * ((1 - cos(2 * pi * x_coordinate(x, y_substrate - 1, x_steps, x_length))) ** M0) \
                                    + p[y_substrate - 3][x]  # plank pg 179 eq 65
            if p[y_substrate - 1][x] - p_old > tolerance or p[y_substrate - 1][x] - p_old < -tolerance:
                in_tol = 0

        # X = Max
        p_old = p[y_substrate - 1][x_steps - 2]
        p[y_substrate - 1][x_steps - 2] = p[y_substrate - 1][x_steps - 3]  # minus 3 because row is even
        if p[y_substrate - 1][x_steps - 2] - p_old > tolerance or p[y_substrate - 1][x_steps - 2] - p_old < -tolerance:
            in_tol = 0

        # Y = Max - 1

        # X = 0
        p_old = p[y_substrate - 2][0]
        p[y_substrate - 2][0] = p[y_substrate - 2][1]  # plank pg 179 eq 70
        if p[y_substrate - 2][0] - p_old > tolerance or p[y_substrate - 2][0] - p_old < -tolerance:
            in_tol = 0

        # Body
        for x in range(1, x_steps - 1, 1):
            p_old = p[y_substrate - 2][x]
            p[y_substrate - 2][x] = K35 * h \
                                    * ((1 - cos(2 * pi * x_coordinate(x, y_substrate - 2, x_steps, x_length))) ** M0) \
                                    + p[y_substrate - 4][x]  # -4 because substrate points are at half mesh points
            if p[y_substrate - 2][x] - p_old > tolerance or p[y_substrate - 2][x] - p_old < -tolerance:
                in_tol = 0  # plank pg 179 eq 65

        # X = Max
        p_old = p[y_substrate - 2][x_steps - 1]
        p[y_substrate - 2][x_steps - 1] = p[y_substrate - 2][x_steps - 2]  # plank pg 179 eq 70
        if p[y_substrate - 2][x_steps - 1] - p_old > tolerance or p[y_substrate - 2][x_steps - 1] - p_old < -tolerance:
            in_tol = 0

        # Interior nodes
        for y in range(y_substrate - 3, 2, -1):

            # X = 0
            p_old = p[y][0]
            p[y][0] = p[y][1]  # plank pg 179 eq 70
            if p[y][0] - p_old > tolerance or p[y][0] - p_old < -tolerance:
                in_tol = 0

            # Body
            if y % 2 == 0:  # If row is even number of substrate mesh points in x is x_steps-1
                for x in range(1, x_steps - 2, 1):  # minus 2 because last mesh point has a boundary condition
                    density = density_scale * ((y_substrate / 2) - 1) \
                              * (occupied_old[y // 2][x] + occupied_old[y // 2][x + 1]) / 2  # Av density right and left
                    p_old = p[y][x]
                    p[y][x] = RELAX1 / (h * h + 2 * K21 * k) \
                              * (0.5 * K21 * k * (p[y][x + 1] + p[y][x - 1] + p[y + 2][x] + p[y - 2][x]
                                                  + pedf[y][x + 1] + pedf[y][x - 1] + pedf[y + 2][x] + pedf[y - 2][x])
                                 + (h * h - 2 * K21 * k - h * h * k * K1 * density / (1 + pedf[y][x])) * pedf[y][x]) \
                              + (1 - RELAX1) * p[y][x]  # plank pg 178 eq 53 crank-nicolson approximation
                    if p[y][x] - p_old > tolerance or p[y][x] - p_old < -tolerance:
                        in_tol = 0

                # X = Max
                p_old = p[y][x_steps - 2]
                p[y][x_steps - 2] = p[y][x_steps - 3]  # plank pg 179 eq 70
                if p[y][x_steps - 2] - p_old > tolerance or p[y][x_steps - 2] - p_old < -tolerance:
                    in_tol = 0

            # Body
            else:
                for x in range(1, x_steps - 1, 1):  # If row is odd number of substrate mesh points in x is nn
                    density = density_scale * ((y_substrate / 2) - 1) \
                              * (occupied_old[(y - 1) // 2][x] + occupied_old[(y + 1) // 2][
                        x]) / 2  # Ave den above and below
                    p_old = p[y][x]
                    p[y][x] = RELAX1 / (h * h + 2 * K21 * k) \
                              * (0.5 * K21 * k * (p[y][x + 1] + p[y][x - 1] + p[y + 2][x] + p[y - 2][x]
                                                  + pedf[y][x + 1] + pedf[y][x - 1] + pedf[y + 2][x] + pedf[y - 2][x])
                                 + (h * h - 2 * K21 * k - h * h * k * K1 * density / (1 + pedf[y][x])) * pedf[y][x]) \
                              + (1 - RELAX1) * p[y][x]
                    if p[y][x] - p_old > tolerance or p[y][x] - p_old < -tolerance:
                        in_tol = 0

                # X = Max
                p_old = p[y][x_steps - 1]
                p[y][x_steps - 1] = p[y][x_steps - 2]  # plank pg 179 eq 70
                if p[y][x_steps - 1] - p_old > tolerance or p[y][x_steps - 1] - p_old < -tolerance:
                    in_tol = 0

        # Y = 2

        # X = 0
        p_old = p[2][0]
        p[2][0] = p[2][1]  # plank pg 179 eq 70
        if p[2][0] - p_old > tolerance or p[2][0] - p_old < -tolerance:
            in_tol = 0

        # Body
        for x in range(1, x_steps - 2, 1):
            p_old = p[2][x]
            p[2][x] = 1 / (K33 + h) * (K33 * p[4][x] + h * pedf_old[0][x])  # plank pg 179, 137 eq 61
            if p[2][x] - p_old > tolerance or p[2][x] - p_old < -tolerance:
                in_tol = 0

        # X = Max
        p_old = p[2][x_steps - 2]
        p[2][x_steps - 2] = p[2][x_steps - 3]  # plank pg 179 eq 70
        if p[2][x_steps - 2] - p_old > tolerance or p[2][x_steps - 2] - p_old < -tolerance:
            in_tol = 0

        # Y = 1: capillary wall

        # X = 0
        p_old = p[1][0]
        p[1][0] = p[1][1]
        if p[1][0] - p_old > tolerance or p[1][0] - p_old < -tolerance:
            in_tol = 0

        # Body
        for x in range(1, x_steps - 1, 1):
            p_old = p[1][x]
            p[1][x] = 1 / (K33 + h) * (K33 * p[3][x] + h * (pedf_old[0][x - 1] + pedf_old[0][x]) / 2)  # Ave
            if p[1][x] - p_old > tolerance or p[1][x] - p_old < -tolerance:
                in_tol = 0

        # X = Max
        p_old = p[1][x_steps - 1]
        p[1][x_steps - 1] = p[1][x_steps - 2]  # plank pg 179 eq 70
        if p[1][x_steps - 1] - p_old > tolerance or p[1][x_steps - 1] - p_old < -tolerance:
            in_tol = 0

    # Cycle through pedf mesh points and set pedf at time step j+1
    for y in range(1, y_substrate, 1):  # If y is even: number of substrate mesh points in x is x_steps - 1
        if y % 2 == 0:
            for x in range(x_steps - 1):
                pedf_old[y][x] = pedf[y][x] * 0
                pedf[y][x] = p[y][x] * 0
                if pedf[y][x] < 0:
                    pedf[y][x] = 0
        else:  # If y is odd: number of substrate mesh points in x is x_steps
            for x in range(x_steps):
                pedf_old[y][x] = pedf[y][x] * 0
                pedf[y][x] = p[y][x] * 0
                if pedf[y][x] < 0:
                    pedf[y][x] = 0

    return pedf, pedf_old
