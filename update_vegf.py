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

    # Initialize v, the substrate matrix used to iterate
    v = zeros((y_substrate, x_steps))
    v_old = zeros((y_substrate, x_steps))
    for y in range(1, y_substrate):
        for x in range(x_steps):
            v[y][x] = vegf[y][x]
            vegf_old[y][x] = vegf[y][x]

    # Capillary
    for x in range(x_steps - 1):
        density = density_scale * (occupied_old[0][x] + occupied_old[0][x+1]) / 2  # Average density
        wall_vegf = (vegf[1][x] + vegf[1][x+1]) / 2  # Average VEGF concentration on capillary wall at time j
        vegf_difference = wall_vegf - vegf[0][x]
        if vegf_difference < 0:
            vegf_difference = 0
        vegf[0][x] = vegf[0][x] + k \
                     * (-K1 * vegf[0][x] * density / (1 + vegf[0][x]) + K2 * vegf_difference)  # plank pg 150 eq 46

    # ECM iteration loop
    in_tol = 0  # tolerance key for the while loop
    while in_tol == 0:
        in_tol = 1  # Set the key to free the loop unless it gets flipped

        # Reset v_old
        for y in range(1, y_substrate):
            for x in range(x_steps):
                v_old[y][x] = v[y][x]

        # Source
        v[y_substrate - 1][0] = v[y_substrate - 1][1]  # plank pg 179 eq 70, x=0
        for x in range(1, x_steps - 2):
            v[y_substrate - 1][x] = K35 * h \
                                    * ((1 - cos(2 * pi * x_coordinate(x, y_substrate - 1, x_steps, x_length))) ** M0) \
                                    + v[y_substrate - 3][x]  # plank pg 179 eq 65, body
        v[y_substrate-1][x_steps-2] = v[y_substrate-1][x_steps-3]  # plank pg 179 eq 70, x=max

        # Source wall (Y = Max - 1)
        v[y_substrate-2][0] = v[y_substrate-2][1]  # plank pg 179 eq 70, x=0
        for x in range(1, x_steps - 1):
            v[y_substrate - 2][x] = K35 * h \
                                    * ((1 - cos(2 * pi * x_coordinate(x, y_substrate - 2, x_steps, x_length))) ** M0) \
                                    + v[y_substrate - 4][x]  # plank pg 179 eq 65, x=body
        v[y_substrate - 2][x_steps - 1] = v[y_substrate - 2][x_steps - 2]  # plank pg 179 eq 70, x=max

        # Interior nodes
        for y in range(y_substrate - 3, 2, -1):
            v[y][0] = v[y][1]  # plank pg 179 eq 70, x=0
            if y % 2 == 0:  # Even rows, body
                for x in range(1, x_steps - 2, 1):
                    density = density_scale * ((y_substrate / 2) - 1) \
                              * (occupied_old[y//2][x] + occupied_old[y//2][x+1]) / 2  # Average density
                    v[y][x] = RELAX1 / (h * h + 2 * K21 * k) \
                              * (0.5 * K21 * k * (v[y][x+1] + v[y][x-1] + v[y+2][x] + v[y-2][x]
                                                  + vegf[y][x+1] + vegf[y][x-1] + vegf[y+2][x] + vegf[y-2][x])
                                 + (h * h - 2 * K21 * k - h * h * k * K1 * density / (1 + vegf[y][x])) * vegf[y][x]) \
                              + (1 - RELAX1) * v[y][x]  # plank pg 178 eq 53
                v[y][x_steps - 2] = v[y][x_steps - 3]  # plank pg 179 eq 70, x=max
            else:  # Odd rows, body
                for x in range(1, x_steps - 1):
                    density = density_scale * ((y_substrate / 2) - 1) \
                              * (occupied_old[(y-1)//2][x] + occupied_old[(y+1)//2][x]) / 2  # Average density
                    v[y][x] = RELAX1 / (h * h + 2 * K21 * k) \
                              * (0.5 * K21 * k * (v[y][x+1] + v[y][x-1] + v[y+2][x] + v[y-2][x]
                                                  + vegf[y][x+1] + vegf[y][x-1] + vegf[y+2][x] + vegf[y-2][x])
                                 + (h * h - 2 * K21 * k - h * h * k * K1 * density / (1 + vegf[y][x])) * vegf[y][x]) \
                              + (1 - RELAX1) * v[y][x]
                v[y][x_steps - 1] = v[y][x_steps - 2]  # plank pg 179 eq 70, x=max

        # Capillary wall (Y = 2, 1)
        v[2][0] = v[2][1]  # plank pg 179 eq 70, x=0
        for x in range(1, x_steps - 2, 1):
            v[2][x] = 1 / (K33 + h) * (K33 * v[4][x] + h * vegf_old[0][x])  # plank pg 179, 137 eq 61, x=body
        v[2][x_steps - 2] = v[2][x_steps - 3]  # plank pg 179 eq 70, x=max
        v[1][0] = v[1][1]  # plank pg 179 eq 70, x=0
        for x in range(1, x_steps - 1):
            v[1][x] = 1 / (K33 + h) * (K33 * v[3][x] + h * (vegf_old[0][x - 1] + vegf_old[0][x]) / 2)  # x=body
        v[1][x_steps - 1] = v[1][x_steps - 2]  # plank pg 179 eq 70, x=max

        # Tolerance Check
        for y in range(1, y_substrate):
            for x in range(x_steps):
                if v[y][x] - v_old[y][x] > tolerance or v[y][x] - v_old[y][x] < -tolerance:
                    in_tol = 0

    # Set vegf at time step j + 1
    for y in range(1, y_substrate, 1):
        for x in range(x_steps):
            vegf[y][x] = v[y][x]

    # VEGF can't go below 0
    for y in range(y_substrate):
        for x in range(x_steps):
            if vegf[y][x] < 0:
                vegf[y][x] = 0

    return
