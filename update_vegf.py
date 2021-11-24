# Description
'''
The purpose of UpdateVEGF is to update the vegf substrate for each new time step
'''


# Imports
# Library Imports
from numpy import zeros
from math import cos
from math import pi

# File Imports
from x_coordinate import x_coordinate

# Values Imports
from the_vault import K1, K2, K35, K21, M0, K33, RELAX1

# Function
def update_vegf(ySubstrate, xSteps, densityScale, occupiedOld, vegf, vegfOld, k, tolerance, h, xLength):

    # Substrate matrix used to iterate
    v = zeros((ySubstrate, xSteps))

    # Update VEGF concentration in capillary
    for x in range(xSteps - 1):

        # Scaled by densityScale to get density instead of the number of cells
        # Average density at cell mesh points to the right and left of the substrate mesh point
        # Use occupiedOld because occupied has already been updated this time step
        density = densityScale * (occupiedOld[0][x] + occupiedOld[0][x + 1]) / 2

        # Average VEGF concentration on capillary wall at time j
        wallVEGFatJ = (vegf[1][x] + vegf[1][x + 1]) / 2
        VEGFdifference = wallVEGFatJ - vegf[0][x]
        if VEGFdifference < 0:
            VEGFdifference = 0
        vegfOld[0][x] = vegf[0][x]

        # Update VEGF concentration in capillary using equation 46
        vegf[0][x] = vegf[0][x] + k * (-K1 * vegf[0][x] * density / (1 + vegf[0][x]) + K2 * VEGFdifference)
        if vegf[0][x] < 0:
            vegf[0][x] = 0

    # Initialize v: value of the previous time step
    for y in range(1, ySubstrate, 1):
        if y % 2 == 0:
            for x in range(xSteps - 1):
                v[y][x] = vegf[y][x]
        else:
            for x in range(xSteps):
                v[y][x] = vegf[y][x]

    # Keep iterating until the value of v at each mesh point changes by less than a certain tolerance
    inTolerance = 0
    while inTolerance == 0:
        inTolerance = 1

    # Update VEGF concentration at boundary at maximum y. Includes source
        # Using equation 70 derivation on page 179 calculate vegf at x = 0
        vOld = v[ySubstrate - 1][0]
        v[ySubstrate - 1][0] = v[ySubstrate - 1][1]

        # If the change in any v and vOld is greater than the tolerance then this loop will continue to run
        if v[ySubstrate - 1][0] - vOld > tolerance or v[ySubstrate - 1][0] - vOld < -tolerance:
            inTolerance = 0

        for x in range(1, xSteps - 2, 1):
            vOld = v[ySubstrate - 1][x]

            # Use EQ 65 and derivation on page 179 to update VEGF concentration at upper boundary
            v[ySubstrate - 1][x] = K35 * h * ((1 - cos(2 * pi * x_coordinate(x, ySubstrate - 1, xSteps, xLength))) ** M0) \
                                   + v[ySubstrate - 3][x]
            if v[ySubstrate - 1][x] - vOld > tolerance or v[ySubstrate - 1][x] - vOld < -tolerance:
                inTolerance = 0

        # Using equation 70 derivation on page 179 calculate vegf at x = max
        vOld = v[ySubstrate - 1][xSteps - 2]
        v[ySubstrate - 1][xSteps - 2] = v[ySubstrate - 1][xSteps - 3]    # minus 3 because row is even
        if v[ySubstrate - 1][xSteps - 2] - vOld > tolerance or v[ySubstrate - 1][xSteps - 2] - vOld < -tolerance:
            inTolerance = 0

    # Update VEGF concentration at boundary at maximum y - 1
        # Using equation 70 derivation on page 179 calculate vegf at x = 0
        vOld = v[ySubstrate - 2][0]
        v[ySubstrate - 2][0] = v[ySubstrate - 2][1]
        if v[ySubstrate - 2][0] - vOld > tolerance or v[ySubstrate - 2][0] - vOld < -tolerance:
            inTolerance = 0

        for x in range(1, xSteps - 1, 1):
            vOld = v[ySubstrate - 2][x]
            # Use EQ 65 and derivation on page 179 to update VEGF concentration at upper boundary - 1
            v[ySubstrate - 2][x] = K35 * h * ((1 - cos(2 * pi * x_coordinate(x, ySubstrate - 2, xSteps, xLength))) ** M0) \
                                   + v[ySubstrate - 4][x]    # -4 because substrate points are at half mesh points
            if v[ySubstrate - 2][x] - vOld > tolerance or v[ySubstrate - 2][x] - vOld < -tolerance:
                inTolerance = 0

        # Using equation 70 derivation on page 179 calculate vegf at x = max
        vOld = v[ySubstrate - 2][xSteps - 1]
        v[ySubstrate - 2][xSteps - 1] = v[ySubstrate - 2][xSteps - 2]
        if v[ySubstrate - 2][xSteps - 1] - vOld > tolerance or v[ySubstrate - 2][xSteps - 1] - vOld < -tolerance:
            inTolerance = 0

    # Cycle through interior rows and calculate new VEGF concentration
        for y in range(ySubstrate - 3, 2, -1):
            # Using equation 70 derivation on page 179 calculate vegf at x = 0
            vOld = v[y][0]
            v[y][0] = v[y][1]
            if v[y][0] - vOld > tolerance or v[y][0] - vOld < -tolerance:
                inTolerance = 0

            # If row is even number of substrate meshpoints in x is xsteps-1
            if y % 2 == 0:
                for x in range(1, xSteps - 2, 1):    # minus 2 because last meshpoint has a boundary condition
                    # DensityScale is squared because density in ECM is per unit area not length
                    # Average density of cell meshpoint to the right and left
                    density = densityScale * ((ySubstrate/2)-1) * (occupiedOld[y // 2][x] + occupiedOld[y // 2][x + 1]) / 2
                    vOld = v[y][x]

                    # Approximate equation 53 using the crank-nicolson method see derivation on page 178
                    v[y][x] = RELAX1 / (h * h + 2 * K21 * k) * (0.5 * K21 * k * (v[y][x + 1] + v[y][x - 1] + v[y + 2][x] + \
                            v[y - 2][x] + vegf[y][x + 1] + vegf[y][x - 1] + vegf[y + 2][x] + vegf[y - 2][x]) + \
                            (h * h - 2 * K21 * k - h * h * k * K1 * density / (1 + vegf[y][x])) * vegf[y][x]) + (1-RELAX1) * v[y][x]
                    if v[y][x] - vOld > tolerance or v[y][x] - vOld < -tolerance:
                        inTolerance = 0

                # Using equation 70 derivation on page 179 calculate vegf at x = max
                vOld = v[y][xSteps - 2]
                v[y][xSteps - 2] = v[y][xSteps - 3]
                if v[y][xSteps - 2] - vOld > tolerance or v[y][xSteps - 2] - vOld < -tolerance:
                    inTolerance = 0

            else:
                # If row is odd number of substrate meshpoints in x is nn
                for x in range(1, xSteps - 1, 1):
                    # Average density of cell meshpoints above and below
                    # y // 2 because there are twice as many points in the y direction because substrate meshpoints are at 1/2
                    density = densityScale * ((ySubstrate/2)-1) * (occupiedOld[(y - 1) // 2][x] + occupiedOld[(y + 1) // 2][x]) / 2
                    vOld = v[y][x]
                    v[y][x] = RELAX1 / (h * h + 2 * K21 * k) * (0.5 * K21 * k * (v[y][x + 1] + v[y][x - 1] + v[y + 2][x] + \
                            v[y - 2][x] + vegf[y][x + 1] + vegf[y][x - 1] + vegf[y + 2][x] + vegf[y - 2][x]) + \
                            (h * h - 2 * K21 * k - h * h * k * K1 * density / (1 + vegf[y][x])) * vegf[y][x]) + (1-RELAX1) * v[y][x]
                    if v[y][x] - vOld > tolerance or v[y][x] - vOld < -tolerance:
                        inTolerance = 0

                # Using equation 70 derivation on page 179 calculate vegf at x = max
                vOld = v[y][xSteps - 1]
                v[y][xSteps - 1] = v[y][xSteps - 2]
                if v[y][xSteps - 1] - vOld > tolerance or v[y][xSteps - 1] - vOld < -tolerance:
                    inTolerance = 0

    # Calculate VEGF concentration at boundary row y = 2
        # Using equation 70 derivation on page 179 calculate vegf at x = 0
        vOld = v[2][0]
        v[2][0] = v[2][1]
        if v[2][0] - vOld > tolerance or v[2][0] - vOld < -tolerance:
            inTolerance = 0

        for x in range(1, xSteps - 2, 1):
            vOld = v[2][x]
            # Use equation 61 see derivation on page 179 of paper and pg 137 of notes
            # Use vegfOld because capillary values have already been updated
            v[2][x] = 1 / (K33 + h) * (K33 * v[4][x] + h * vegfOld[0][x])
            if v[2][x] - vOld > tolerance or v[2][x] - vOld < -tolerance:
                inTolerance = 0

        # Using equation 70 derivation on page 179 calculate vegf at x = max
        vOld = v[2][xSteps - 2]
        v[2][xSteps - 2] = v[2][xSteps - 3]
        if v[2][xSteps - 2] - vOld > tolerance or v[2][xSteps - 2] - vOld < -tolerance:
            inTolerance = 0

    # Calculate VEGF concentratin at boundary row y = 1: capillary wall
        vOld = v[1][0]
        v[1][0] = v[1][1]

        if v[1][0] - vOld > tolerance or v[1][0] - vOld < -tolerance:
            inTolerance = 0
        for x in range(1, xSteps - 1, 1):
            vOld = v[1][x]

            # Take average because there is no mesh point in the capillary directly below
            v[1][x] = 1 / (K33 + h) * (K33 * v[3][x] + h * (vegfOld[0][x - 1] + vegfOld[0][x]) / 2)
            if v[1][x] - vOld > tolerance or v[1][x] - vOld < -tolerance:
                inTolerance = 0

        # Using equation 70 derivation on page 179 calculate vegf at x = max
        vOld = v[1][xSteps - 1]
        v[1][xSteps - 1] = v[1][xSteps - 2]
        if v[1][xSteps - 1] - vOld > tolerance or v[1][xSteps - 1] - vOld < -tolerance:
            inTolerance = 0

    # Cycle through VEGF mesh points and set VEGF at time step j+1
    for y in range(1, ySubstrate, 1):
        # If y is even: number of substrate mesh points in x is xSteps - 1
        if y % 2 == 0:
            for x in range(xSteps - 1):
                vegfOld[y][x] = vegf[y][x]
                vegf[y][x] = v[y][x]
                if vegf[y][x] < 0:
                    vegf[y][x] = 0

        # If y is odd: number of substrate mesh points in x is xSteps
        else:
            for x in range(xSteps):
                vegfOld[y][x] = vegf[y][x]
                vegf[y][x] = v[y][x]
                if vegf[y][x] < 0:
                    vegf[y][x] = 0

    return
