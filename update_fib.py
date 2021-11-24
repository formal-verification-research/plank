# Description
'''
The purpose of UpdateFibronectin is to update the fibronectin substrate for each new time step
'''


# Imports
# Library Imports
from numpy import zeros

# Values Imports
from the_vault import K4, K5, K6, K22, K23, RELAX2


# Function
def update_fib(ySubstrate, xSteps, densityScale, occupiedOld, fibronectin, fibronectinOld, k, protease,
                      tolerance, h, pedfOld):

    # Substrate matrix used to iterate
    f = zeros((ySubstrate, xSteps))

    # Update fibronectin concentration in capillary
    for x in range(xSteps-1):

        # Scaled by densityScale to get density instead of the number of cells
        # Average density at cell mesh points to the right and left of the substrate mesh point
        # Use occupiedOld because occupied has already been updated this time step
        density = densityScale * (occupiedOld[0][x] + occupiedOld[0][x+1]) / 2
        fibronectinOld[0][x] = fibronectin[0][x]

        # Calculate fibronectin concentration in capillary using equation 48 and equation 51
        fibronectin[0][x] = fibronectin[0][x] + k * (K4 * fibronectin[0][x] * (1 - fibronectin[0][x]) * density - K5 * \
                            protease[0][x] * fibronectin[0][x] / (1 + pedfOld[0][x] + K6 * fibronectin[0][x]))

        # Make sure fibronectin never goes negative
        if fibronectin[0][x] < 0:
            fibronectin[0][x] = 0

        # One is the initial and max value of the fibronectin concentration
        if fibronectin[0][x] > 1:
            fibronectin[0][x] = 1

    # Initialize f: value of previous time step
    for y in range(1, ySubstrate, 1):
        if y % 2 == 0:
            for x in range(xSteps-1):
                f[y][x] = fibronectin[y][x]
        else:
            for x in range(xSteps):
                f[y][x] = fibronectin[y][x]

    # Keep iterating until the value of f at each meshpoint changes by less than a certain tolerance
    inTolerance = 0
    while inTolerance == 0:
        inTolerance = 1

        # Update fibronectin concentration at boundary at maximum y. Includes source
        # Using equation 71 derivation on page 179 calculate fib at x = 0
        fOld = f[ySubstrate - 1][0]
        f[ySubstrate - 1][0] = f[ySubstrate - 1][1]

        # If the change in any f and fOld is ever greater than the tolerance then this loop will continue to run
        if f[ySubstrate - 1][0] - fOld > tolerance or f[ySubstrate - 1][0] - fOld < -tolerance:
            inTolerance = 0

        for x in range(1, xSteps-2, 1):
            fOld = f[ySubstrate - 1][x]

            # Use EQ 66 and derivation on page 179 to update fib concentration at upper boundary
            f[ySubstrate - 1][x] = f[ySubstrate - 3][x]
            if f[ySubstrate - 1][x] - fOld > tolerance or f[ySubstrate - 1][x] - fOld < -tolerance:
                inTolerance = 0

        # Using equation 71 derivation on page 179 calculate fib at x = max
        fOld = f[ySubstrate - 1][xSteps - 2]
        f[ySubstrate - 1][xSteps - 2] = f[ySubstrate - 1][xSteps - 3]
        if f[ySubstrate - 1][xSteps - 2] - fOld > tolerance or f[ySubstrate - 1][xSteps - 2] - fOld < -tolerance:
            inTolerance = 0

        # Update fib concentration at boundary at maximum y - 1
        # Using equation 71 derivation on page 179 calculate fib at x = 0
        fOld = f[ySubstrate-2][0]
        f[ySubstrate-2][0] = f[ySubstrate-2][1]
        if f[ySubstrate-2][0] - fOld > tolerance or f[ySubstrate-2][0] - fOld < -tolerance:
            inTolerance = 0

        # Use EQ 66 and derivation on page 179 to update fib concentration at upper boundary - 1
        for x in range(1, xSteps-1, 1):
            fOld = f[ySubstrate-2][x]
            f[ySubstrate-2][x] = f[ySubstrate-4][x]
            if f[ySubstrate-2][x] - fOld > tolerance or f[ySubstrate-2][x] - fOld < -tolerance:
                inTolerance = 0

        # Using equation 71 derivation on page 179 calculate fib at x = max
        fOld = f[ySubstrate-2][xSteps-1]
        f[ySubstrate-2][xSteps-1] = f[ySubstrate-2][xSteps-2]
        if f[ySubstrate-2][xSteps-1] - fOld > tolerance or f[ySubstrate-2][xSteps-1] - fOld < -tolerance:
            inTolerance = 0

        # Cycle through interior rows and calculate new fibronectin concentration
        for y in range(ySubstrate-3, 2, -1):
            # Using equation 71 derivation on page 179 calculate fib at x = 0
            fOld = f[y][0]
            f[y][0] = f[y][1]
            if f[y][0] - fOld > tolerance or f[y][0] - fOld < -tolerance:
                inTolerance = 0

            # If row is even number of substrate meshpoints in x is nn-1
            if y % 2 == 0:
                for x in range(1, xSteps-2, 1):    # minus 2 because last meshpoint has a boundary condition
                    fOld = f[y][x]

                    # Approximate equation 55 and equation 59 using the crank-nicolson method see derivation on page 179
                    f[y][x] = RELAX2 * k / (h*h+2*K22*k) * (0.5*K22 * (f[y][x+1]+f[y][x-1]+f[y+2][x]+f[y-2][x] + \
                            fibronectin[y][x+1]+fibronectin[y][x-1]+fibronectin[y+2][x]+fibronectin[y-2][x]) + \
                            (h*h/k - 2*K22 + K23*h*h*(1-fibronectin[y][x]) - K5*h*h*protease[y][x]/\
                             (1+pedfOld[y][x]+K6*fibronectin[y][x])) * fibronectin[y][x]) + (1-RELAX2) * f[y][x]

                    if f[y][x] - fOld > tolerance or f[y][x] - fOld < -tolerance:
                        inTolerance = 0

                # Using equation 71 derivation on page 179 calculate fib at x = max
                fOld = f[y][xSteps-2]
                f[y][xSteps-2] = f[y][xSteps-3]
                if f[y][xSteps-2] - fOld > tolerance or f[y][xSteps-2] - fOld < -tolerance:
                    inTolerance = 0

            else:
                # If row is odd number of substrate meshpoints in x is nn
                for x in range(1, xSteps-1, 1):
                    fOld = f[y][x]

                    # Approximate equation 55 and equation 59 using the crank-nicolson method see derivation on page 179
                    f[y][x] = RELAX2 * k / (h * h + 2 * K22 * k) * (0.5 * K22 * (f[y][x+1] + f[y][x-1] + f[y+2][x] + f[y-2][x] \
                            + fibronectin[y][x+1] + fibronectin[y][x-1] + fibronectin[y+2][x] + fibronectin[y-2][x]) \
                            + (h * h / k - 2 * K22 + K23 * h * h * (1 - fibronectin[y][x]) - K5 * h * h * \
                            protease[y][x] / (1 + pedfOld[y][x] + K6 * fibronectin[y][x]))* fibronectin[y][x]) + (1-RELAX2) * f[y][x]

                    if f[y][x] - fOld > tolerance or f[y][x] - fOld < -tolerance:
                        inTolerance = 0

                # Using equation 71 derivation on page 179 calculate fib at x = max
                fOld = f[y][xSteps-1]
                f[y][xSteps-1] = f[y][xSteps-2]
                if f[y][xSteps-1] - fOld > tolerance or f[y][xSteps-1] - fOld < -tolerance:
                    inTolerance = 0

        # Calculate fibronectin concentration at boundary row y = 2
        # Using equation 71 derivation on page 179 calculate fib at x = 0
        fOld = f[2][0]
        f[2][0] = f[2][1]
        if f[2][0] - fOld > tolerance or f[2][0] - fOld < -tolerance:
            inTolerance=0

        # Use EQ 62 and derivation on page 179 to update fib concentration at lower boundary
        for x in range(1, xSteps-2, 1):
            fOld = f[2][x]
            f[2][x] = f[4][x]
            if f[2][x] - fOld > tolerance or f[2][x] - fOld < -tolerance:
                inTolerance = 0

        # Using equation 71 derivation on page 179 calculate fib at x = max
        fOld = f[2][xSteps-2]
        f[2][xSteps-2] = f[2][xSteps-3]
        if f[2][xSteps-2] - fOld > tolerance or f[2][xSteps-2] - fOld < -tolerance:
            inTolerance = 0

        # Calculate fibronectin concentration at boundary row y = 1: capillary wall
        # Using equation 71 derivation on page 179 calculate fib at x = 0
        fOld = f[1][0]
        f[1][0] = f[1][1]
        if f[1][0] - fOld > tolerance or f[1][0] - fOld < -tolerance:
            inTolerance=0

        # Use EQ 62 and derivation on page 179 to update fib concentration at lower boundary
        for x in range(1, xSteps-1, 1):
            fOld = f[1][x]
            f[1][x] = f[3][x]
            if f[1][x] - fOld > tolerance or f[1][x] - fOld < -tolerance:
                inTolerance = 0

        # Using equation 71 derivation on page 179 calculate fib at x = max
        fOld = f[1][xSteps-1]
        f[1][xSteps-1] = f[1][xSteps-2]
        if f[1][xSteps-1] - fOld > tolerance or f[1][xSteps-1] - fOld < -tolerance:
            inTolerance = 0

    # Cycle through substrate meshpoints and set fibronectin at time step j+1
    for y in range(1, ySubstrate, 1):
        if y % 2 == 0:
            for x in range(xSteps-1):
                fibronectinOld[y][x] = fibronectin[y][x]
                fibronectin[y][x] = f[y][x]
                if fibronectin[y][x] < 0:
                    fibronectin[y][x] = 0
                if fibronectin[y][x] > 1:
                    fibronectin[y][x] = 1
        else:
            for x in range(xSteps):
                fibronectinOld[y][x] = fibronectin[y][x]
                fibronectin[y][x] = f[y][x]
                if fibronectin[y][x] < 0:
                    fibronectin[y][x] = 0
                if fibronectin[y][x] > 1:
                    fibronectin[y][x] = 1

    return
