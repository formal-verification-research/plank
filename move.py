# Description
# move physically moves the EC to its location for the current time step


# Function
def move(cell, time, stay, left, right, up, random_num, y_position, x_position, occupied):
    # Don't need to update x_position because I assumed in the beginning that
    # x_position[cell][time+1] = x_position[cell][time]

    # Down
    if random_num > stay + left + right + up:
        y_position[cell][time+1] = y_position[cell][time] + 1

    # Up
    elif random_num > stay + left + right:
        y_position[cell][time+1] = y_position[cell][time] - 1

    # Right
    elif random_num > stay + left:
        x_position[cell][time+1] = x_position[cell][time] + 1

    # Left
    elif random_num > stay:
        x_position[cell][time+1] = x_position[cell][time] - 1

    # Update occupancy based on cell movement
    if random_num > stay:
        occupied[y_position[cell][time]][x_position[cell][time]] -= 1
        occupied[y_position[cell][time+1]][x_position[cell][time+1]] += 1
