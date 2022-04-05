# Description
# tau calculates movement probability based off of substrate values


# Imports
from parameter_vault import K6, K10, K14, K27, K31, K11, K15, K28, K32, K12, K29, K13, K30
from parameter_vault import GAMMA1, GAMMA2, GAMMA3, GAMMA4, GAMMA5, GAMMA6


# Function
def tau(c, f, v, y):

    ACTIVEC = c / (1 + K6 * f)  # Active protease

    if y == 0:  # plank eq 50
        return (((ACTIVEC + K10) / (ACTIVEC + K11)) ** GAMMA1) * (((f + K12) / (f + K13)) ** GAMMA2) \
               * (((v + K14) / (v + K15)) ** GAMMA3)
    else:  # plank eq 58
        return (((ACTIVEC + K27) / (ACTIVEC + K28)) ** GAMMA4) * (((f + K29) / (f + K30)) ** GAMMA5) \
               * (((v + K31) / (v + K32)) ** GAMMA6)
