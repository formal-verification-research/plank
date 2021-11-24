# Description
# x_coordinate gives the dimensionless value of the distance along the x direction


# Function
def x_coordinate(x, y, x_steps, x_length):
    if y % 2 == 0:
        return ((x + 0.5) / x_steps) * x_length
    else:
        return (x / x_steps) * x_length
