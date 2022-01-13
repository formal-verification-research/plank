# Description
# update_fib updates the fib substrate for each new time step


# Imports
from numpy import zeros
from the_vault import K4, K5, K6, K22, K23, RELAX2


# Function
def update_fib(y_substrate, x_steps, density_scale, occupied_old, fib, fib_old, k, pro, tolerance, h):

    # Initialize f, the matrix used for iteration
    f = zeros((y_substrate, x_steps))
    f_old = zeros((y_substrate, x_steps))
    for y in range(y_substrate):
        for x in range(x_steps):
            f[y][x] = fib[y][x]
            fib_old[y][x] = fib[y][x]

    # Capillary
    for x in range(x_steps - 1):
        density = density_scale * (occupied_old[0][x] + occupied_old[0][x+1]) / 2  # Average density
        ca = pro[0][x] / (1 + K6 * fib[0][x])  # plank pg 151 eq 51
        fib[0][x] = fib[0][x] + k * ((K4 * fib[0][x]
                                      * (1 - fib[0][x]) * density) - (K5 * ca * fib[0][x]))  # plank pg 150 eq 48

    # ECM iteration loop
    in_tol = 0  # Tolerance key for the while loop
    while in_tol == 0:
        in_tol = 1  # Set the key to free the loop unless it gets flipped

        # Reset f_old
        for y in range(1, y_substrate):
            for x in range(x_steps):
                f_old[y][x] = f[y][x]

        # Source
        f[y_substrate-1][0] = f[y_substrate-1][1]  # plank pg 179 eq 71, x=0
        for x in range(1, x_steps - 2):
            f[y_substrate-1][x] = f[y_substrate-3][x]  # plank pg 179 eq 66, body
        f[y_substrate - 1][x_steps - 2] = f[y_substrate - 1][x_steps - 3]  # plank pg 179 eq 71, x=max

        # Source wall (Y = max - 1)
        f[y_substrate-2][0] = f[y_substrate-2][1]  # plank pg 179 eq 71, x=0
        for x in range(1, x_steps - 1):
            f[y_substrate-2][x] = f[y_substrate-4][x]  # plank pg 179 eq 66, x=body
        f[y_substrate-2][x_steps-1] = f[y_substrate-2][x_steps-2]  # plank pg 179 eq 71, x=max

        # Interior nodes
        for y in range(y_substrate - 3, 2, -1):
            f[y][0] = f[y][1]  # plank pg 179 eq 71, x=0
            if y % 2 == 0:  # Even rows, body
                for x in range(1, x_steps - 2, 1):
                    f[y][x] = RELAX2 * k / (h * h + 2 * K22 * k) \
                              * (0.5 * K22 * (f[y][x+1] + f[y][x-1] + f[y+2][x] + f[y-2][x]
                                              + fib[y][x+1] + fib[y][x-1] + fib[y+2][x] + fib[y-2][x])
                                 + (h * h / k - 2 * K22 + K23 * h * h * (1 - fib[y][x])
                                    - K5 * h * h * pro[y][x] / (1 + K6 * fib[y][x])) * fib[y][x]) \
                              + (1 - RELAX2) * f[y][x]  # plank pg 179 eq 55, 59
                f[y][x_steps-2] = f[y][x_steps-3]  # plank pg 179 eq 71, x=max
            else:  # Odd rows, body
                for x in range(1, x_steps - 1):
                    f[y][x] = RELAX2 * k / (h * h + 2 * K22 * k) \
                              * (0.5 * K22 * (f[y][x+1] + f[y][x-1] + f[y+2][x] + f[y-2][x]
                                              + fib[y][x+1] + fib[y][x-1] + fib[y+2][x] + fib[y-2][x])
                                 + (h * h / k - 2 * K22 + K23 * h * h * (1 - fib[y][x])
                                    - K5 * h * h * pro[y][x] / (1 + K6 * fib[y][x])) * fib[y][x]) \
                              + (1-RELAX2) * f[y][x]  # plank pg 179 eq 55, 59
                f[y][x_steps-1] = f[y][x_steps-2]  # plank pg 179 eq 71, x=max

        # Capillary wall (Y = 2, 1)
        f[2][0] = f[2][1]  # plank pg 179 eq 71, x=0
        for x in range(1, x_steps - 2, 1):
            f[2][x] = f[4][x]  # plank pg 179 eq 62, x=body
        f[2][x_steps-2] = f[2][x_steps-3]  # plank pg 179 eq 71, x=max
        f[1][0] = f[1][1]  # plank pg 179 eq 71, x=0
        for x in range(1, x_steps - 1, 1):
            f[1][x] = f[3][x]  # plank pg 179 eq 62, x=body
        f[1][x_steps-1] = f[1][x_steps-2]  # plank pg 179 eq 71, x=max

        # Tolerance Check
        for y in range(1, y_substrate):
            for x in range(x_steps):
                if f[y][x] - f_old[y][x] > tolerance or f[y][x] - f_old[y][x] < -tolerance:
                    in_tol = 0

    # Set fib at time step j + 1
    for y in range(1, y_substrate, 1):
        for x in range(x_steps):
            fib[y][x] = f[y][x]

    # Fib cannot rise above 1 or go below 0
    for y in range(y_substrate):
        for x in range(x_steps):
            if fib[y][x] > 1:
                fib[y][x] = 1
            if fib[y][x] < 0:
                fib[y][x] = 0

    return
