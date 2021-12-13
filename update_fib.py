# Description
# update_fib updates the fib substrate for each new time step


# Imports
from numpy import zeros
from the_vault import K4, K5, K6, K22, K23, RELAX2


# Function
def update_fib(y_substrate, x_steps, density_scale, occupied_old, fib, fib_old, k, pro, tolerance, h):

    # Initialize f
    f = zeros((y_substrate, x_steps))
    for y in range(1, y_substrate, 1):
        if y % 2 == 0:
            for x in range(x_steps-1):
                f[y][x] = fib[y][x]
        else:
            for x in range(x_steps):
                f[y][x] = fib[y][x]

    # Capillary
    for x in range(x_steps - 1):
        density = density_scale * (occupied_old[0][x] + occupied_old[0][x+1]) / 2  # Ave density at right and left nodes 
        fib_old[0][x] = fib[0][x]
        fib[0][x] = fib[0][x] + k * (K4 * fib[0][x] * (1 - fib[0][x]) * density - K5 * pro[0][x]  # plank paper eq 48
                                     * fib[0][x] / (1 + K6 * fib[0][x]))  # plank paper eq 51
        if fib[0][x] < 0:  
            fib[0][x] = 0  # fib doesn't drop below 0
        if fib[0][x] > 1:
            fib[0][x] = 1  # fib doesn't rise above 1

    # ECM
    in_tol = 0  # tolerance key for the while loop
    while in_tol == 0:
        in_tol = 1

        # Y = Max

        # X = 0
        f_old = f[y_substrate - 1][0]
        f[y_substrate - 1][0] = f[y_substrate - 1][1]  # plank pg 179 eq 71
        if f[y_substrate - 1][0] - f_old > tolerance or f[y_substrate - 1][0] - f_old < -tolerance:  # tolerance check
            in_tol = 0

        # Body
        for x in range(1, x_steps-2, 1):
            f_old = f[y_substrate - 1][x]
            f[y_substrate - 1][x] = f[y_substrate - 3][x]  # plank pg 179 eq 66
            if f[y_substrate - 1][x] - f_old > tolerance or f[y_substrate - 1][x] - f_old < -tolerance:  # tolerance check
                in_tol = 0

        # X = Max
        f_old = f[y_substrate - 1][x_steps - 2]
        f[y_substrate - 1][x_steps - 2] = f[y_substrate - 1][x_steps - 3]  # plank pg 179 eq 71
        if f[y_substrate - 1][x_steps - 2] - f_old > tolerance or f[y_substrate - 1][x_steps - 2] - f_old < -tolerance:
            in_tol = 0  # tolerance check

        # Y = Max - 1

        # X = 0
        f_old = f[y_substrate-2][0]
        f[y_substrate-2][0] = f[y_substrate-2][1]  # plank pg 179 eq 71
        if f[y_substrate-2][0] - f_old > tolerance or f[y_substrate-2][0] - f_old < -tolerance:
            in_tol = 0

        # Body
        for x in range(1, x_steps-1, 1):
            f_old = f[y_substrate-2][x]
            f[y_substrate-2][x] = f[y_substrate-4][x]  # plank pg 179 eq 66
            if f[y_substrate-2][x] - f_old > tolerance or f[y_substrate-2][x] - f_old < -tolerance:
                in_tol = 0

        # X = Max
        f_old = f[y_substrate-2][x_steps-1]
        f[y_substrate-2][x_steps-1] = f[y_substrate-2][x_steps-2]  # plank pg 179 eq 71
        if f[y_substrate-2][x_steps-1] - f_old > tolerance or f[y_substrate-2][x_steps-1] - f_old < -tolerance:
            in_tol = 0

        # Interior Nodes
        for y in range(y_substrate-3, 2, -1):

            # X = 0
            f_old = f[y][0]
            f[y][0] = f[y][1]  # plank pg 179 eq 71
            if f[y][0] - f_old > tolerance or f[y][0] - f_old < -tolerance:
                in_tol = 0

            # Body
            if y % 2 == 0:  # If row is even number of substrate nodes in x is nn-1
                for x in range(1, x_steps-2, 1):
                    f_old = f[y][x]
                    f[y][x] = RELAX2 * k / (h * h + 2 * K22 * k) \
                              * (0.5 * K22 * (f[y][x+1] + f[y][x-1] + f[y+2][x] + f[y-2][x]
                                              + fib[y][x+1] + fib[y][x-1] + fib[y+2][x] + fib[y-2][x])
                                 + (h * h / k - 2 * K22 + K23 * h * h * (1 - fib[y][x]) - K5 * h * h * pro[y][x]
                                    / (1 + K6 * fib[y][x])) * fib[y][x]) + (1 - RELAX2) * f[y][x]  # plank pg 179 eq 55, 59
                    if f[y][x] - f_old > tolerance or f[y][x] - f_old < -tolerance:
                        in_tol = 0

                # X = Max
                f_old = f[y][x_steps-2]
                f[y][x_steps-2] = f[y][x_steps-3]  # plank pg 179 eq 71
                if f[y][x_steps-2] - f_old > tolerance or f[y][x_steps-2] - f_old < -tolerance:
                    in_tol = 0

            # Body
            else:  # If row is odd number of substrate nodes in x is nn
                for x in range(1, x_steps-1, 1):
                    f_old = f[y][x]
                    f[y][x] = RELAX2 * k / (h * h + 2 * K22 * k) \
                              * (0.5 * K22 * (f[y][x+1] + f[y][x-1] + f[y+2][x] + f[y-2][x]
                                              + fib[y][x+1] + fib[y][x-1] + fib[y+2][x] + fib[y-2][x])
                                 + (h * h / k - 2 * K22 + K23 * h * h * (1 - fib[y][x]) - K5 * h * h * pro[y][x]
                                    / (1 + K6 * fib[y][x])) * fib[y][x]) + (1-RELAX2) * f[y][x]  # plank pg 179 eq 55, 59
                    if f[y][x] - f_old > tolerance or f[y][x] - f_old < -tolerance:
                        in_tol = 0

                # X = Max
                f_old = f[y][x_steps-1]
                f[y][x_steps-1] = f[y][x_steps-2]  # plank pg 179 eq 71
                if f[y][x_steps-1] - f_old > tolerance or f[y][x_steps-1] - f_old < -tolerance:
                    in_tol = 0

        # Y = 2

        # X = 0
        f_old = f[2][0]
        f[2][0] = f[2][1]  # plank pg 179 eq 71
        if f[2][0] - f_old > tolerance or f[2][0] - f_old < -tolerance:
            in_tol = 0

        # Body
        for x in range(1, x_steps-2, 1):
            f_old = f[2][x]
            f[2][x] = f[4][x]  # plank pg 179 eq 62
            if f[2][x] - f_old > tolerance or f[2][x] - f_old < -tolerance:
                in_tol = 0

        # X = Max
        f_old = f[2][x_steps-2]
        f[2][x_steps-2] = f[2][x_steps-3]  # plank pg 179 eq 71
        if f[2][x_steps-2] - f_old > tolerance or f[2][x_steps-2] - f_old < -tolerance:
            in_tol = 0

        # Y = 1

        # X = 0
        f_old = f[1][0]
        f[1][0] = f[1][1]  # plank pg 179 eq 71
        if f[1][0] - f_old > tolerance or f[1][0] - f_old < -tolerance:
            in_tol = 0

        # Body
        for x in range(1, x_steps-1, 1):
            f_old = f[1][x]
            f[1][x] = f[3][x]  # plank pg 179 eq 62
            if f[1][x] - f_old > tolerance or f[1][x] - f_old < -tolerance:
                in_tol = 0

        # X = Max
        f_old = f[1][x_steps-1]
        f[1][x_steps-1] = f[1][x_steps-2]  # plank pg 179 eq 71
        if f[1][x_steps-1] - f_old > tolerance or f[1][x_steps-1] - f_old < -tolerance:
            in_tol = 0

    # Cycle through substrate mesh points and set fib at time step j+1
    for y in range(1, y_substrate, 1):
        if y % 2 == 0:
            for x in range(x_steps-1):
                fib_old[y][x] = fib[y][x]
                fib[y][x] = f[y][x]
                if fib[y][x] < 0:
                    fib[y][x] = 0
                if fib[y][x] > 1:
                    fib[y][x] = 1
        else:
            for x in range(x_steps):
                fib_old[y][x] = fib[y][x]
                fib[y][x] = f[y][x]
                if fib[y][x] < 0:
                    fib[y][x] = 0
                if fib[y][x] > 1:
                    fib[y][x] = 1

    return
