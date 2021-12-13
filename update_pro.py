# Description
# update_pro updates the pro substrate for each new time step


# Imports
from the_vault import K1, K3


# Function
def update_pro(y_substrate, x_steps, density_scale, occupied_old, pro, pro_old, k, vegf_old, pedf_old):

    # Capillary
    for x in range(x_steps - 1):
        density = density_scale * (occupied_old[0][x] + occupied_old[0][x+1]) / 2  # Ave density at right and left nodes
        pro_old[0][x] = pro[0][x]
        pro[0][x] = pro[0][x] + k * (K1 * (vegf_old[0][x] - pedf_old[0][x])  # plank eq 47
                                     * density / (vegf_old[0][x] + 1) - K3 * pro[0][x])
        if pro[0][x] < 0:
            pro[0][x] = 0

    # ECM

    # Y = 1
    for x in range(x_steps):
        density = density_scale * (y_substrate / 2 - 1) * occupied_old[1][x]
        pro_old[1][x] = pro[1][x]
        pro[1][x] = pro[1][x] + k * (K1 * (vegf_old[1][x] - pedf_old[1][x])  # plank eq 57
                                     * density / (vegf_old[1][x] + 1) - K3 * pro[1][x])
        if pro[1][x] < 0:
            pro[1][x] = 0

    # Interior Nodes
    for y in range(2, y_substrate, 1):

        if y % 2 == 0:  # If y is even, number of substrate mesh points in x-direction is x_steps-1
            for x in range(x_steps - 1):
                density = density_scale * (y_substrate / 2 - 1) \
                          * (occupied_old[y // 2][x] + occupied_old[y // 2][x + 1]) / 2  # Ave density right and left
                pro_old[y][x] = pro[y][x]
                pro[y][x] = pro[y][x] + k * (K1 * (vegf_old[y][x] - pedf_old[y][x])  # plank eq 54
                                             * density / (vegf_old[y][x] + 1) - K3 * pro[y][x])
                if pro[y][x] < 0:
                    pro[y][x] = 0

        else:  # If y is odd, number of substrate mesh points in x-direction is x_steps
            for x in range(x_steps):
                # Average density at cell mesh points above and below substrate mesh point at time j
                density = density_scale * (y_substrate / 2 - 1) \
                          * (occupied_old[(y-1) // 2][x] + occupied_old[(y+1) // 2][x]) / 2  # Av density right and left
                pro_old[y][x] = pro[y][x]
                pro[y][x] = pro[y][x] + k * (K1 * (vegf_old[y][x] - pedf_old[y][x])
                                             * density / (vegf_old[y][x] + 1) - K3 * pro[y][x])  # plank eq 54
                if pro[y][x] < 0:
                    pro[y][x] = 0

    return
