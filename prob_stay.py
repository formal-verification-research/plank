"""
Description:
The file prob_stay determines the probability that the EC will not move for the current time step
"""


# Function
def prob_stay(y, lam, k):

    # The equation that dictates the probability of staying is found on Plank Paper Pages 152-153
    if y == 0:
        stay = 1 - 2 * lam * k
    else:
        stay = 1 - 4 * lam * k

    return stay
