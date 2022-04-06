"""
Description:
The tau file calculates EC movement probability in a single direction based off of the substrate values
"""


# Imports
from parameter_vault import K6, K10, K14, K27, K31, K11, K15, K28, K32, K12, K29, K13, K30
from parameter_vault import GAMMA1, GAMMA2, GAMMA3, GAMMA4, GAMMA5, GAMMA6


# Function
def tau(c, f, v, y):

    # Only the active Protease counts for the movement probabilities
    ACTIVEC = c / (1 + K6 * f)

    # Probability calculation based on Protease, Fibronectin, and VEGF (Plank Paper Equations 50 and 58)
    if y == 0:
        return (((ACTIVEC + K10) / (ACTIVEC + K11)) ** GAMMA1) * (((f + K12) / (f + K13)) ** GAMMA2) \
               * (((v + K14) / (v + K15)) ** GAMMA3)
    else:
        return (((ACTIVEC + K27) / (ACTIVEC + K28)) ** GAMMA4) * (((f + K29) / (f + K30)) ** GAMMA5) \
               * (((v + K31) / (v + K32)) ** GAMMA6)
