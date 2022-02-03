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

    # Initialize p, the substrate matrix used to iterate
    p = zeros((y_substrate, x_steps))
    p_old = zeros((y_substrate, x_steps))
    for y in range(1, y_substrate, 1):
        for x in range(x_steps):
            p[y][x] = pedf[y][x]
            pedf_old[y][x] = pedf[y][x]

    # Capillary
    for x in range(x_steps - 1):
        density = density_scale * (occupied_old[0][x] + occupied_old[0][x + 1]) / 2  # Average density
        wall_pedf = (pedf[1][x] + pedf[1][x+1]) / 2  # Average pedf concentration on capillary wall at time j
        pedf_difference = wall_pedf - pedf[0][x]
        if pedf_difference < 0:
            pedf_difference = 0
        pedf[0][x] = pedf[0][x] + k \
                     * (-K1 * pedf[0][x] * density / (1 + pedf[0][x]) + K2 * pedf_difference)  # plank pg 150 eq 46

    # ECM iteration loop
    in_tol = 0  # tolerance key for the while loop
    while in_tol == 0:
        in_tol = 1  # Set the key to free the loop unless it gets flipped

        # Reset p_old
        for y in range(1, y_substrate):
            for x in range(x_steps):
                p_old[y][x] = p[y][x]

        # Source
        p[y_substrate - 1][0] = p[y_substrate - 1][1]  # plank pg 179 eq 70, x=0
        for x in range(1, x_steps - 2):
            p[y_substrate - 1][x] = K35 * h \
                                    * ((1 - cos(2 * pi * x_coordinate(x, y_substrate - 1, x_steps, x_length))) ** M0) \
                                    + p[y_substrate - 3][x]  # plank pg 179 eq 65, body
        p[y_substrate - 1][x_steps - 2] = p[y_substrate - 1][x_steps - 3]  # plank pg 179 eq 70, x=max

        # Source wall (Y = Max - 1)
        p[y_substrate - 2][0] = p[y_substrate - 2][1]  # plank pg 179 eq 70, x=0
        for x in range(1, x_steps - 1):
            p[y_substrate - 2][x] = K35 * h \
                                    * ((1 - cos(2 * pi * x_coordinate(x, y_substrate - 2, x_steps, x_length))) ** M0) \
                                    + p[y_substrate - 4][x]  # plank pg 179 eq 65, x=body
        p[y_substrate - 2][x_steps - 1] = p[y_substrate - 2][x_steps - 2]  # plank pg 179 eq 70, x=max

        # Interior nodes
        for y in range(y_substrate - 3, 2, -1):
            p[y][0] = p[y][1]  # plank pg 179 eq 70, x=0
            if y % 2 == 0:  # Even rows, body
                for x in range(1, x_steps - 2, 1):
                    density = density_scale * ((y_substrate / 2) - 1) \
                              * (occupied_old[y//2][x] + occupied_old[y//2][x+1]) / 2  # Average density
                    p[y][x] = RELAX1 / (h * h + 2 * K21 * k) \
                              * (0.5 * K21 * k * (p[y][x+1] + p[y][x-1] + p[y+2][x] + p[y-2][x]
                                                  + pedf[y][x+1] + pedf[y][x-1] + pedf[y+2][x] + pedf[y-2][x])
                                 + (h * h - 2 * K21 * k - h * h * k * K1 * density / (1 + pedf[y][x])) * pedf[y][x]) \
                              + (1 - RELAX1) * p[y][x]  # plank pg 178 eq 53
                p[y][x_steps - 2] = p[y][x_steps - 3]  # plank pg 179 eq 70, x=max
            else:  # Odd rows, body
                for x in range(1, x_steps - 1):
                    density = density_scale * ((y_substrate / 2) - 1) \
                              * (occupied_old[(y-1)//2][x] + occupied_old[(y+1)//2][x]) / 2  # Average density
                    p[y][x] = RELAX1 / (h * h + 2 * K21 * k) \
                              * (0.5 * K21 * k * (p[y][x+1] + p[y][x-1] + p[y+2][x] + p[y-2][x]
                                                  + pedf[y][x+1] + pedf[y][x-1] + pedf[y+2][x] + pedf[y-2][x])
                                 + (h * h - 2 * K21 * k - h * h * k * K1 * density / (1 + pedf[y][x])) * pedf[y][x]) \
                              + (1 - RELAX1) * p[y][x]
                p[y][x_steps - 1] = p[y][x_steps - 2]  # plank pg 179 eq 70, x=max

        # Capillary wall (Y = 2, 1)
        p[2][0] = p[2][1]  # plank pg 179 eq 70, x=0
        for x in range(1, x_steps - 2, 1):
            p[2][x] = 1 / (K33 + h) * (K33 * p[4][x] + h * pedf_old[0][x])  # plank pg 179, 137 eq 61, x=body
        p[2][x_steps - 2] = p[2][x_steps - 3]  # plank pg 179 eq 70, x=max
        p[1][0] = p[1][1]  # plank pg 179 eq 70, x=0
        for x in range(1, x_steps - 1):
            p[1][x] = 1 / (K33 + h) * (K33 * p[3][x] + h * (pedf_old[0][x - 1] + pedf_old[0][x]) / 2)  # x=body
        p[1][x_steps - 1] = p[1][x_steps - 2]  # plank pg 179 eq 70, x=max

        # Tolerance Check
        for y in range(1, y_substrate):
            for x in range(x_steps):
                if p[y][x] - p_old[y][x] > tolerance or p[y][x] - p_old[y][x] < -tolerance:
                    in_tol = 0

    # Set pedf at time step j + 1
    for y in range(1, y_substrate, 1):
        for x in range(x_steps):
            pedf[y][x] = p[y][x] * 0

    # PEDF can't go below 0
    for y in range(y_substrate):
        for x in range(x_steps):
            if pedf[y][x] < 0:
                pedf[y][x] = 0

    return pedf, pedf_old
