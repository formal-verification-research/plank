# Description
'''
The purpose of pMove is to calculate the odds of the EC traveling to each space around it
'''


# Imports
# File Import
from Tau import tau


# Function
def pMove(x, y, direction, protease, fibronectin, vegf, pedf, xSteps, ySteps, lam, k):

    # T is the chance that it will go left, right, up, and down
    T = [0, 0, 0, 0]

    # Left
    if x > 0:
        T[0] = tau(protease[2*y][x-1], fibronectin[2*y][x-1], vegf[2*y][x-1], pedf[2*y][x-1], y)
    else:
        T[0] = 0

    # Right
    # Don't need to add 1 to x because the x substrate should be offset by a half step from the cell matrix meshpoints
    if x < xSteps-1:
        T[1] = tau(protease[2*y][x], fibronectin[2*y][x], vegf[2*y][x], pedf[2*y][x], y)
    else:
        T[1] = 0

    # Up
    if y > 1:
        T[2] = tau(protease[2*y-1][x], fibronectin[2*y-1][x], vegf[2*y-1][x], pedf[2*y-1][x], y)
    else:
        T[2] = 0

    # Down
    if 0 < y < ySteps - 1:
        T[3] = tau(protease[2*y+1][x], fibronectin[2*y+1][x], vegf[2*y+1][x], pedf[2*y+1][x], y)
    else:
        T[3] = 0

    # Probability equations from page 153 and 152 of the plank paper
    if y == 0:
        pmove = 2 * lam * k * T[direction] / (T[0]+T[1])
    else:
        pmove = 4 * lam * k * T[direction] / (T[0]+T[1]+T[2]+T[3])

    return pmove, T
