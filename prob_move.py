# Description
# prob_move calculates the odds of the EC traveling to each space around it


# Imports
from tau import tau


# Function
def prob_move(x, y, direction, pro, fib, vegf, pedf, x_steps, y_steps, lam, k):
    # T is the chance that it will go left, right, up, and down
    T = [0, 0, 0, 0]

    # Left
    if x > 0:
        T[0] = tau(pro[2*y][x-1], fib[2*y][x-1], vegf[2*y][x-1], pedf[2*y][x-1], y)
    else:
        T[0] = 0

    # Right, don't need to add 1 because the x substrate is offset a half step from the cell matrix mesh points
    if x < x_steps-1:
        T[1] = tau(pro[2*y][x], fib[2*y][x], vegf[2*y][x], pedf[2*y][x], y)
    else:
        T[1] = 0

    # Up
    if y > 1:
        T[2] = tau(pro[2*y-1][x], fib[2*y-1][x], vegf[2*y-1][x], pedf[2*y-1][x], y)
    else:
        T[2] = 0

    # Down
    if 0 < y < y_steps - 1:
        T[3] = tau(pro[2*y+1][x], fib[2*y+1][x], vegf[2*y+1][x], pedf[2*y+1][x], y)
    else:
        T[3] = 0

    if y == 0:
        pmove = 2 * lam * k * T[direction] / (T[0]+T[1])  # plank paper page 152 and 153
    else:
        pmove = 4 * lam * k * T[direction] / (T[0]+T[1]+T[2]+T[3])

    return pmove, T
