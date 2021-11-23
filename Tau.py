# Description
'''
The purpose of Tau is to calculate movement probability based off of substrate values
'''


# Imports
# Values Imports
from the_vault import K6, K10, K14, K27, K31, K11, K15, K28, K32, K12, K29, K13, K30
from the_vault import GAMMA1, GAMMA2, GAMMA3, GAMMA4, GAMMA5, GAMMA6


# Function
def tau(c, f, v, p, y):

    # This is the active protease
    activeC = c / (1 + p + K6 * f)

    if y == 0:
        # equation 50
        return (((activeC + K10) / (activeC + K11)) ** GAMMA1) * (((f + K12) / (f + K13)) ** GAMMA2) * (
                    ((v + K14) / (v + K15)) ** GAMMA3)
    else:
        # equation 58
        return (((activeC + K27) / (activeC + K28)) ** GAMMA4) * (((f + K29) / (f + K30)) ** GAMMA5) * (
                    ((v + K31) / (v + K32)) ** GAMMA6)
