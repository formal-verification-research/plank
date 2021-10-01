# Description
'''
The purpose of pStay is to determine the chances that the cell will not more for the current time step
'''


# Function
def pStay(y, lam, k):

    # Probability equations found on page 152 and 153 of the plank paper
    if y == 0:
        pstay = 1 - 2 * lam * k
    else:
        pstay = 1 - 4 * lam * k

    return pstay
