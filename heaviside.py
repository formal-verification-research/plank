# Description
# heaviside returns a 1 or a 0 based on the input, usually cell maturity or substrate timings


# Function
def heaviside(variable):
    if variable >= 0:
        return 1
    else:
        return 0
