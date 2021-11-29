# Description
'''
The purpose of Updatepedf is to update the pedf substrate for each new time step
'''


# Imports
from numpy import zeros
from math import cos
from math import pi
from the_vault import V0
from x_coordinate import x_coordinate
from the_vault import K1, K2, K35, K21, M0, K33, RELAX1
from heaviside import heaviside


# Function
def update_pedf(ySubstrate, xSteps, densityScale, occupiedOld, pedf, pedfOld, k, tolerance, h, xLength,
                current_time_step):

    # Substrate matrix used to iterate
    v = zeros((ySubstrate, xSteps))

    # Update pedf concentration in capillary
    for x in range(xSteps - 1):

        # Scaled by densityScale to get density instead of the number of cells
        # Average density at cell mesh points to the right and left of the substrate mesh point
        # Use occupiedOld because occupied has already been updated this time step
        density = densityScale * (occupiedOld[0][x] + occupiedOld[0][x + 1]) / 2

        # Average pedf concentration on capillary wall at time j
        wallpedfatJ = (pedf[1][x] + pedf[1][x + 1]) / 2
        pedfdifference = wallpedfatJ - pedf[0][x]
        if pedfdifference < 0:
            pedfdifference = 0
        pedfOld[0][x] = pedf[0][x]

        # Update pedf concentration in capillary using equation 46
        pedf[0][x] = pedf[0][x] + k * (-K1 * pedf[0][x] * density / (1 + pedf[0][x]) + K2 * pedfdifference)
        if pedf[0][x] < 0:
            pedf[0][x] = 0

    # Initialize v: value of the previous time step
    for y in range(1, ySubstrate, 1):
        if y % 2 == 0:
            for x in range(xSteps - 1):
                v[y][x] = pedf[y][x]
        else:
            for x in range(xSteps):
                v[y][x] = pedf[y][x]

    # Keep iterating until the value of v at each mesh point changes by less than a certain tolerance
    inTolerance = 0
    while inTolerance == 0:
        inTolerance = 1

    # Update pedf concentration at boundary at maximum y. Includes source
        # Using equation 70 derivation on page 179 calculate pedf at x = 0
        vOld = v[ySubstrate - 1][0]
        v[ySubstrate - 1][0] = v[ySubstrate - 1][1]

        # If the change in any v and vOld is greater than the tolerance then this loop will continue to run
        if v[ySubstrate - 1][0] - vOld > tolerance or v[ySubstrate - 1][0] - vOld < -tolerance:
            inTolerance = 0

        for x in range(1, xSteps - 2, 1):
            vOld = v[ySubstrate - 1][x]

            # Use EQ 65 and derivation on page 179 to update pedf concentration at upper boundary
            v[ySubstrate - 1][x] = 1
            if v[ySubstrate - 1][x] - vOld > tolerance or v[ySubstrate - 1][x] - vOld < -tolerance:
                inTolerance = 0

        # Using equation 70 derivation on page 179 calculate pedf at x = max
        vOld = v[ySubstrate - 1][xSteps - 2]
        v[ySubstrate - 1][xSteps - 2] = v[ySubstrate - 1][xSteps - 3]    # minus 3 because row is even
        if v[ySubstrate - 1][xSteps - 2] - vOld > tolerance or v[ySubstrate - 1][xSteps - 2] - vOld < -tolerance:
            inTolerance = 0

    # Update pedf concentration at boundary at maximum y - 1
        # Using equation 70 derivation on page 179 calculate pedf at x = 0
        vOld = v[ySubstrate - 2][0]
        v[ySubstrate - 2][0] = v[ySubstrate - 2][1]
        if v[ySubstrate - 2][0] - vOld > tolerance or v[ySubstrate - 2][0] - vOld < -tolerance:
            inTolerance = 0

        for x in range(1, xSteps - 1, 1):
            vOld = v[ySubstrate - 2][x]
            # Use EQ 65 and derivation on page 179 to update pedf concentration at upper boundary - 1
            v[ySubstrate - 2][x] = 1  # -4 because substrate points are at half mesh points
            if v[ySubstrate - 2][x] - vOld > tolerance or v[ySubstrate - 2][x] - vOld < -tolerance:
                inTolerance = 0

        # Using equation 70 derivation on page 179 calculate pedf at x = max
        vOld = v[ySubstrate - 2][xSteps - 1]
        v[ySubstrate - 2][xSteps - 1] = v[ySubstrate - 2][xSteps - 2]
        if v[ySubstrate - 2][xSteps - 1] - vOld > tolerance or v[ySubstrate - 2][xSteps - 1] - vOld < -tolerance:
            inTolerance = 0

    # Cycle through interior rows and calculate new pedf concentration
        for y in range(ySubstrate - 3, 2, -1):
            # Using equation 70 derivation on page 179 calculate pedf at x = 0
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
                            v[y - 2][x] + pedf[y][x + 1] + pedf[y][x - 1] + pedf[y + 2][x] + pedf[y - 2][x]) + \
                            (h * h - 2 * K21 * k - h * h * k * K1 * density / (1 + pedf[y][x])) * pedf[y][x]) + (1-RELAX1) * v[y][x]
                    if v[y][x] - vOld > tolerance or v[y][x] - vOld < -tolerance:
                        inTolerance = 0

                # Using equation 70 derivation on page 179 calculate pedf at x = max
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
                            v[y - 2][x] + pedf[y][x + 1] + pedf[y][x - 1] + pedf[y + 2][x] + pedf[y - 2][x]) + \
                            (h * h - 2 * K21 * k - h * h * k * K1 * density / (1 + pedf[y][x])) * pedf[y][x]) + (1-RELAX1) * v[y][x]
                    if v[y][x] - vOld > tolerance or v[y][x] - vOld < -tolerance:
                        inTolerance = 0

                # Using equation 70 derivation on page 179 calculate pedf at x = max
                vOld = v[y][xSteps - 1]
                v[y][xSteps - 1] = v[y][xSteps - 2]
                if v[y][xSteps - 1] - vOld > tolerance or v[y][xSteps - 1] - vOld < -tolerance:
                    inTolerance = 0

    # Calculate pedf concentration at boundary row y = 2
        # Using equation 70 derivation on page 179 calculate pedf at x = 0
        vOld = v[2][0]
        v[2][0] = v[2][1]
        if v[2][0] - vOld > tolerance or v[2][0] - vOld < -tolerance:
            inTolerance = 0

        for x in range(1, xSteps - 2, 1):
            vOld = v[2][x]
            # Use equation 61 see derivation on page 179 of paper and pg 137 of notes
            # Use pedfOld because capillary values have already been updated
            v[2][x] = 1 / (K33 + h) * (K33 * v[4][x] + h * pedfOld[0][x])
            if v[2][x] - vOld > tolerance or v[2][x] - vOld < -tolerance:
                inTolerance = 0

        # Using equation 70 derivation on page 179 calculate pedf at x = max
        vOld = v[2][xSteps - 2]
        v[2][xSteps - 2] = v[2][xSteps - 3]
        if v[2][xSteps - 2] - vOld > tolerance or v[2][xSteps - 2] - vOld < -tolerance:
            inTolerance = 0

    # Calculate pedf concentratin at boundary row y = 1: capillary wall
        vOld = v[1][0]
        v[1][0] = v[1][1]

        if v[1][0] - vOld > tolerance or v[1][0] - vOld < -tolerance:
            inTolerance = 0
        for x in range(1, xSteps - 1, 1):
            vOld = v[1][x]

            # Take average because there is no mesh point in the capillary directly below
            v[1][x] = 1 / (K33 + h) * (K33 * v[3][x] + h * (pedfOld[0][x - 1] + pedfOld[0][x]) / 2)
            if v[1][x] - vOld > tolerance or v[1][x] - vOld < -tolerance:
                inTolerance = 0

        # Using equation 70 derivation on page 179 calculate pedf at x = max
        vOld = v[1][xSteps - 1]
        v[1][xSteps - 1] = v[1][xSteps - 2]
        if v[1][xSteps - 1] - vOld > tolerance or v[1][xSteps - 1] - vOld < -tolerance:
            inTolerance = 0

    # Cycle through pedf mesh points and set pedf at time step j+1
    for y in range(1, ySubstrate, 1):
        # If y is even: number of substrate mesh points in x is xSteps - 1
        if y % 2 == 0:
            for x in range(xSteps - 1):
                pedfOld[y][x] = pedf[y][x]
                pedf[y][x] = v[y][x]
                if pedf[y][x] < 0:
                    pedf[y][x] = 0

        # If y is odd: number of substrate mesh points in x is xSteps
        else:
            for x in range(xSteps):
                pedfOld[y][x] = pedf[y][x]
                pedf[y][x] = v[y][x]
                if pedf[y][x] < 0:
                    pedf[y][x] = 0

    return
