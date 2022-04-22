"""
Description
The update_fib file updates the Fibronectin substrate for each new time step
"""


# Imports
from numpy import zeros
from parameter_vault import K4, K5, K6, K22, K23, RELAX2, x_steps, tolerance


# Function
def update_fib(y_substrate, density_cap, ec_old, fib, fib_old, k, pro, h):

    # Update the Fibronectin concentration inside the parent blood vessel using Plank Eq 48, 51 Pgs 150-151
    for x in range(x_steps-1):
        fib_old[0][x] = fib[0][x]
        density = density_cap * (ec_old[0][x] + ec_old[0][x+1]) / 2
        ca = pro[0][x] / (1 + K6 * fib[0][x])
        fib[0][x] = fib_old[0][x] + k \
                    * ((K4 * fib_old[0][x] * (1 - fib_old[0][x]) * density) - (K5 * ca * fib_old[0][x]))

    # Set Fibronectin for the previous time step
    for y in range(1, y_substrate):
        for x in range(x_steps):
            fib_old[y][x] = fib[y][x]

    # Create a Fibronectin array to help with the ECM equations
    f = zeros((y_substrate, x_steps))
    for y in range(1, y_substrate):
        for x in range(x_steps):
            f[y][x] = fib[y][x]

    # Update the Fibronectin in the capillary wall and in the ECM using an iterative loop
    tol = 0
    while tol == 0:
        tol = 1

        # Take the current iteration and place it into the fib array
        for y in range(1, y_substrate):
            for x in range(x_steps):
                fib[y][x] = f[y][x]

        # Update the Fibronectin at the RPE layer and surrounding row using Plank Eq 66, 71, Pg 179
        for x in range(1, x_steps-2):
            f[y_substrate-1][x] = f[y_substrate-3][x]
        f[y_substrate-1][0] = f[y_substrate-1][1]
        f[y_substrate-1][x_steps-2] = f[y_substrate-1][x_steps-3]
        for x in range(1, x_steps-1):
            f[y_substrate-2][x] = f[y_substrate-4][x]
        f[y_substrate-2][0] = f[y_substrate-2][1]
        f[y_substrate-2][x_steps-1] = f[y_substrate-2][x_steps-2]

        # Update the Fibronectin at the interior nodes of the ECM using Plank Eq 55, 59, 71, Pg 179
        for y in range(y_substrate-3, 2, -1):
            if y % 2 == 0:
                for x in range(1, x_steps-2):
                    f[y][x] = RELAX2 * k / (h * h + 2 * K22 * k) \
                              * (K22 * 0.5 * (f[y][x+1] + f[y][x-1] + f[y+2][x] + f[y-2][x]
                                            + fib_old[y][x+1] + fib_old[y][x-1] + fib_old[y+2][x] + fib_old[y-2][x])
                                 + (h * h / k - 2 * K22 + K23 * h * h * (1 - fib_old[y][x])
                                    - K5 * h * h * pro[y][x] / (1 + K6 * fib_old[y][x]) / (1 + K6 * fib_old[y][x]))
                                 * fib_old[y][x]) + (1 - RELAX2) * f[y][x]
                f[y][0] = f[y][1]
                f[y][x_steps-2] = f[y][x_steps-3]
            else:
                for x in range(1, x_steps-1):
                    f[y][x] = RELAX2 * k / (h * h + 2 * K22 * k) \
                              * (K22 * 0.5 * (f[y][x+1] + f[y][x-1] + f[y+2][x] + f[y-2][x]
                                            + fib_old[y][x+1] + fib_old[y][x-1] + fib_old[y+2][x] + fib_old[y-2][x])
                                 + (h * h / k - 2 * K22 + K23 * h * h * (1 - fib_old[y][x])
                                    - K5 * h * h * pro[y][x] / (1 + K6 * fib_old[y][x]) / (1 + K6 * fib_old[y][x]))
                                 * fib_old[y][x]) + (1-RELAX2) * f[y][x]
                f[y][0] = f[y][1]
                f[y][x_steps-1] = f[y][x_steps-2]

        # Update the Fibronectin inside the blood vessel wall using Plank Eq 62, 71 Pg 179
        for x in range(1, x_steps-2):
            f[2][x] = f[4][x]
        f[2][0] = f[2][1]
        f[2][x_steps-2] = f[2][x_steps-3]
        for x in range(1, x_steps-1):
            f[1][x] = f[3][x]
        f[1][0] = f[1][1]
        f[1][x_steps-1] = f[1][x_steps-2]

        # Check to make sure each point meets the tolerance before the while loop breaks
        for y in range(1, y_substrate):
            for x in range(x_steps):
                if f[y][x] - fib[y][x] > tolerance or f[y][x] - fib[y][x] < -tolerance:
                    tol = 0

    # Cycle through Fibronectin mesh points and set Fibronectin for the future time step
    for y in range(1, y_substrate):
        for x in range(x_steps):
            fib[y][x] = f[y][x]

    # Fib cannot rise above 1 or go below 0
    for y in range(y_substrate):
        for x in range(x_steps):
            if fib[y][x] > 1:
                fib[y][x] = 1
            if fib[y][x] < 0:
                fib[y][x] = 0

    return fib, fib_old
