"""
Description:
The x_coordinate file returns the non-dimensionalized value of the distance along the x axis
"""


# Function
def x_coordinate(x, y, x_steps, x_length):

    # Distance along the x axis (0-1)
    if y % 2 == 0:
        return ((x + 0.5) / x_steps) * x_length
    else:
        return (x / x_steps) * x_length
