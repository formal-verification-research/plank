"""
Description
The update_pro file updates the Protease substrate for each new time step. This file also incorporates PEDF into the
model by assuming that the EC no longer emits Protease based on VEGF alone. Now the EC samples the surrounding VEGF and
the PEDF and finds the difference, and uses that value to secrete Protease. This assumes that PEDF desensitizes the EC
to VEGF and makes it less able to create Protease based on surrounding VEGF.
"""


# Imports
from parameter_vault import K1, K3, x_steps


# Function
def update_pro(y_substrate, density_cap, density_ecm, ec_old, pro, pro_old, k, vegf_old, pedf_old, current_time_step, total_number_time_steps):

    # Place the Protease into the previous time step array
    for y in range(y_substrate):
        for x in range(x_steps):
            pro_old[y][x] = pro[y][x]

    # Update the Protease in the parent blood vessel using Plank Eq 47 Pg 150
    for x in range(x_steps-1):
        density = density_cap * (ec_old[0][x] + ec_old[0][x+1]) / 2
        desensitize = vegf_old[0][x] - pedf_old[0][x]
        if desensitize < 0:
            desensitize = 0
        pro[0][x] = pro_old[0][x] + k * (K1 * desensitize * density / (vegf_old[0][x] + 1) - K3 * pro_old[0][x])

    # Update the Protease in the parent blood vessel wall using Plank Eq 54 Pg 151
    for x in range(x_steps):
        density = density_ecm * ec_old[1][x]
        desensitize = vegf_old[1][x] - pedf_old[1][x]

        if desensitize < 0:
            desensitize = 0
        pro[1][x] = pro_old[1][x] + k * (K1 * desensitize * density / (vegf_old[1][x] + 1) - K3 * pro_old[1][x])

    # Update the Protease in the interior nodes in the ECM using Plank Eq 54 Pg 151
    for y in range(2, y_substrate):
        if y % 2 == 0:
            for x in range(x_steps-1):
                density = density_ecm * (ec_old[y//2][x] + ec_old[y//2][x+1]) / 2
                desensitize = vegf_old[y][x] - pedf_old[y][x]
                if desensitize < 0:
                    desensitize = 0
                pro[y][x] = pro_old[y][x] + k * (K1 * desensitize * density / (vegf_old[y][x] + 1) - K3 * pro_old[y][x])
        else:
            for x in range(x_steps):
                density = density_ecm * (ec_old[(y-1)//2][x] + ec_old[(y+1)//2][x]) / 2
                desensitize = vegf_old[y][x] - pedf_old[y][x]
                if desensitize < 0:
                    desensitize = 0
                pro[y][x] = pro_old[y][x] + k * (K1 * desensitize * density / (vegf_old[y][x] + 1) - K3 * pro_old[y][x])

    # Make sure the Protease cannot be a negative value
    for y in range(y_substrate):
        for x in range(x_steps):
            if pro[y][x] < 0:
                pro[y][x] = 0

    return pro, pro_old
