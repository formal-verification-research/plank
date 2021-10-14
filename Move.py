# Description
'''
The purpose of Move is to physically move the EC to the location for the next time step
'''


# Function
def move(cell, time, stay, left, right, up, randomNumber, yPosition, xPosition, occupied):

    # Don't need to update xPosition because I assumed in the beginning that
    # xPosition[cell][time+1] = xPosition[cell][time]
    # Down
    if randomNumber > stay + left + right + up:
        yPosition[cell][time+1] = yPosition[cell][time] + 1

    # Up
    elif randomNumber > stay + left + right:
        yPosition[cell][time+1] = yPosition[cell][time] - 1

    # Right
    elif randomNumber > stay + left:
        xPosition[cell][time+1] = xPosition[cell][time] + 1

    # Left
    elif randomNumber > stay:
        xPosition[cell][time+1] = xPosition[cell][time] - 1

    # Update occupancy based on cell movement
    if randomNumber > stay:
        occupied[yPosition[cell][time]][xPosition[cell][time]] -= 1
        occupied[yPosition[cell][time+1]][xPosition[cell][time+1]] += 1
