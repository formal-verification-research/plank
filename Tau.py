# Description
'''
The purpose of Tau is to calculate movement probability based off of substrate values
'''


# Imports
# Values Imports
from Values import K6, K10, K14, K27, K31, K11, K15, K28, K32, K12, K29, K13, K30
from Values import gamma1, gamma2, gamma3, gamma4, gamma5, gamma6


# Function
def tau(c, f, v, p, y):

    # This is the active protease
    activeC = c / (1 + K6 * f)

    if y == 0:
        # equation 50
        return (((activeC + K10) / (activeC + K11)) ** gamma1) * (((f + K12) / (f + K13)) ** gamma2) * (
                    (((v-p) + K14) / ((v-p) + K15)) ** gamma3)
    else:
        # equation 58
        return (((activeC + K27) / (activeC + K28)) ** gamma4) * (((f + K29) / (f + K30)) ** gamma5) * (
                    (((v-p) + K31) / ((v-p) + K32)) ** gamma6)
