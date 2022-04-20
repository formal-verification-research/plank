# Description
# update_pro updates the pro substrate for each new time step


# Imports
from parameter_vault import K1, K3


# Function
def update_pro(y_substrate, x_steps, density_scale, occupied_old, pro, pro_old, k, vegf_old, pedf_old):

    # Set the previous time step
    for y in range(y_substrate):
        for x in range(x_steps):
            pro_old[y][x] = pro[y][x]

    # Capillary
    for x in range(x_steps - 1):
        density = density_scale * (occupied_old[0][x] + occupied_old[0][x+1]) / 2  # Average density
        pro[0][x] = pro_old[0][x] + k * (K1 * (vegf_old[0][x] - pedf_old[0][x]) * density / (vegf_old[0][x] + 1)
                                         - K3 * pro_old[0][x])  # plank pg 150 eq 47

    # Capillary Wall
    for x in range(x_steps):
        density = density_scale * (y_substrate / 2 - 1) * occupied_old[1][x]  # Average density
        pro[1][x] = pro_old[1][x] + k * (K1 * (vegf_old[1][x] - pedf_old[1][x]) * density / (vegf_old[1][x] + 1)
                                         - K3 * pro_old[1][x])  # plank pg 151 eq 54

    # Interior Nodes
    for y in range(2, y_substrate, 1):

        # Even rows
        if y % 2 == 0:
            for x in range(x_steps - 1):
                density = density_scale * (y_substrate / 2 - 1) * (occupied_old[y//2][x]
                                                                   + occupied_old[y//2][x+1]) / 2  # Average density
                pro[y][x] = pro_old[y][x] + k \
                            * (K1 * (vegf_old[y][x] - pedf_old[y][x]) * density / (vegf_old[y][x] + 1)
                               - K3 * pro_old[y][x])  # plank pg 151 eq 54

        # Odd rows
        else:
            for x in range(x_steps):
                density = density_scale * (y_substrate / 2 - 1) * (occupied_old[(y-1)//2][x]
                                                                   + occupied_old[(y+1)//2][x]) / 2  # Average density
                pro[y][x] = pro_old[y][x] + k \
                            * (K1 * (vegf_old[y][x] - pedf_old[y][x]) * density / (vegf_old[y][x] + 1)
                               - K3 * pro_old[y][x])  # plank pg 151 eq 54

    # Protease cannot be negative
    for y in range(y_substrate):
        for x in range(x_steps):
            if pro[y][x] < 0:
                pro[y][x] = 0

    return pro, pro_old
