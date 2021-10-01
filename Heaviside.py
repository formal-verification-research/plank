# Description
'''
The purpose of Heaviside is to return an all-or-nothing value
'''


# Function
def heaviside(variable):

    if variable >= 0:
        return 1
    else:
        return 0
