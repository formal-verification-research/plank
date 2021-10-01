# Description
'''
The purpose of Lamda is to return a ration of the matrix size to the size between meshpoints
'''


# Function
def lamda(xLength, h):

    # non-dimensionalized equation on page 150
    # can't use lambda because that is a python keyword
    lam = (xLength ** 2) / (h ** 2)

    return lam
