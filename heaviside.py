"""
Description:
The heaviside file is meant to take an input and return an "all-or-nothing" value of 0 or 1. It is used to decide if
enough time has passed before the EC is able to divide again.
"""


# Function
def heaviside(variable):

    # Returning that "all-or-nothing" value
    if variable >= 0:
        return 1
    else:
        return 0
