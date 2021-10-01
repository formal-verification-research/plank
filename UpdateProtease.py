# Description
'''
The purpose of UpdateProtease is to update the protease substrate for each new time step
'''


# Imports
# Values Imports
from Values import K1, K3


# Function
def updateProtease(ySubstrate, xSteps, densityScale, occupiedOld, protease, proteaseOld, k, vegfOld):

    # Update protease concentration in capillary
    for x in range(xSteps -1):

        # Average density at cell meshpoints to the right and left of substrate meshpoint at time j
        density = densityScale * (occupiedOld[0][x] + occupiedOld[0][ x +1]) / 2
        proteaseOld[0][x] = protease[0][x]

        # Use equation 47 to update protease concentration in the capillary
        # Use vegfOld because vegf has already been updated this time step
        protease[0][x] = protease[0][x] + k * (K1 * vegfOld[0][x] * density / (vegfOld[0][x ] +1) - K3 * protease[0][x])
        if protease[0][x] < 0:
            protease[0][x] = 0

    # Update protease concentration at boundary y = 1
    for x in range(xSteps):

        # don't take into account cells still in the capillary
        density = densityScale * (ySubstrate / 2 -1) * occupiedOld[1][x]
        proteaseOld[1][x] = protease[1][x]

        # use equation 54 (same as 47 but for the ECM instead of the capillary)
        protease[1][x] = protease[1][x] + k * (K1 * vegfOld[1][x] * density / (vegfOld[1][x ] +1) - K3 * protease[1][x])
        if protease[1][x] < 0:
            protease[1][x] = 0

    # Cycle through substrate meshpoints
    for y in range(2, ySubstrate, 1):

        # If y is even, number of substrate meshpoints in x-direction is xSteps-1
        if y % 2 == 0:
            for x in range(xSteps -1):

                # Average density at cell meshpoints to the right and left of substrate meshpoint at time j
                density = densityScale * (ySubstrate / 2 -1) * (occupiedOld[y // 2][x] + occupiedOld[y // 2][x + 1]) / 2
                proteaseOld[y][x] = protease[y][x]

                # Use equation 54 (same as 47 but for the ECM instead of the capillary)
                protease[y][x] = protease[y][x] + k * (K1 * vegfOld[y][x] * density / (vegfOld[y][x] + 1) - K3 * protease[y][x])
                if protease[y][x] < 0:
                    protease[y][x] = 0

        # If y is odd, number of substrate meshpoints in x-direction is xSteps
        else:
            for x in range(xSteps):
                # Average density at cell meshpoints above and below substrate meshpoint at time j
                density = densityScale * (ySubstrate / 2 - 1) * (occupiedOld[(y - 1) // 2][x] + occupiedOld[(y + 1) // 2][x]) / 2
                proteaseOld[y][x] = protease[y][x]

                # use equation 54 (same as 47 but for the ECM instead of the capillary)
                protease[y][x] = protease[y][x] + k * (K1 * vegfOld[y][x] * density / (vegfOld[y][x] + 1) - K3 * protease[y][x])
                if protease[y][x] < 0:
                    protease[y][x] = 0

    return
