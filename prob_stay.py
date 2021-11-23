# Description
# prob_stay determines the chance that the cell will not move for the current time step


# Function
def prob_stay(y, lam, k):

    if y == 0:
        pstay = 1 - 2 * lam * k  # Probability equations found on page 152 and 153 of the plank paper
    else:
        pstay = 1 - 4 * lam * k

    return pstay
