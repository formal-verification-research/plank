"""
Description:
The move function physically moves the EC to its new location for the current time step
"""

# Function
def move(cell, time, stay, left, right, up, random_num, y_position, x_position, ec):

    # You only need to update the position of the EC on the axis that it travelled, because we have already assumed the
    # EC stayed put for this current time step, so the other axis that remains unchanged is already set.

    # Move the EC in the direction that the number generator picked (Down, Up, Right, Left is the order coded here)
    if random_num > stay + left + right + up:
        y_position[cell][time+1] = y_position[cell][time] + 1
    elif random_num > stay + left + right:
        y_position[cell][time+1] = y_position[cell][time] - 1
    elif random_num > stay + left:
        x_position[cell][time+1] = x_position[cell][time] + 1
    elif random_num > stay:
        x_position[cell][time+1] = x_position[cell][time] - 1

    # Update the EC array based on the cell movement for this time step
    if random_num > stay:
        ec[y_position[cell][time]][x_position[cell][time]] -= 1
        ec[y_position[cell][time+1]][x_position[cell][time+1]] += 1

    return y_position, x_position, ec
