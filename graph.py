# Description
# graph turns the different substrate matrices and the workspace matrix into a graph for better visual 
# interpretation of the results


# Imports
from matplotlib import pyplot as plt
from mpl_toolkits import mplot3d
from matplotlib import cm
from the_vault import V1, F0, VE, L, DP


# Function
def graph(y_substrate, x_steps, vegf, pedf, fib, pro, x_vector, y_vector, workspace, current_time_step):

    vegf_z = []  # vegf height data for the 3D graph
    pedf_z = []  # pedf height data for the 3D graph
    fib_z = []  # fib height data for the 3D graph
    pro_z = []  # pro height data for the 3D graph

    # Add the appropriate info for odd and even lines to the 3D vector
    for y in range(y_substrate):
        if y % 2 == 0:
            for x in range(x_steps-1):
                vegf_z.append(vegf[y][x] / V1)
                pedf_z.append(pedf[y][x] / V1)
                fib_z.append(fib[y][x] * F0)
                pro_z.append(pro[y][x] / V1)
        else:
            for x in range(x_steps):
                vegf_z.append(vegf[y][x] / V1)
                pedf_z.append(pedf[y][x] / V1)
                fib_z.append(fib[y][x] * F0)
                pro_z.append(pro[y][x] / V1)

    # Create the overall file
    file_name = "Time = " + str(current_time_step * L * L / DP * 60) + " Minutes" + ".png"
    file = open(file_name, "a+")
    fig = plt.figure(figsize=plt.figaspect(10.0))

    # Create the EC color map
    ax = fig.add_subplot(5, 1, 1)
    ax.imshow(workspace, cmap='Purples')
    ax.title.set_text('Angiogenesis')

    # Create the vegf 3D graph
    ax = fig.add_subplot(5, 1, 2, projection='3d')
    ax.plot_trisurf(x_vector, y_vector, vegf_z, cmap='hot', edgecolor='none')
    # ax.set_zlim3d(0, 10)
    ax.title.set_text('VEGF')
    ax.set_zlabel('uM')

    # Create the pedf 3D graph
    ax = fig.add_subplot(5, 1, 3, projection='3d')
    ax.plot_trisurf(x_vector, y_vector, pedf_z, cmap='cool', edgecolor='none')
    # ax.set_zlim3d(0, 10)
    ax.title.set_text('PEDF')
    ax.set_zlabel('uM')

    # Create the fib 3D graph
    ax = fig.add_subplot(5, 1, 4, projection='3d')
    ax.plot_trisurf(x_vector, y_vector, fib_z, cmap='plasma', edgecolor='none')
    # ax.set_zlim3d(0, 1)
    ax.title.set_text('Fibronectin')
    ax.set_zlabel('uM')

    # Create the protease 3D graph
    ax = fig.add_subplot(5, 1, 5, projection='3d')
    ax.plot_trisurf(x_vector, y_vector, pro_z, cmap='viridis', edgecolor='none')
    # ax.set_zlim3d(0, 10)
    ax.title.set_text('Protease')
    ax.set_zlabel('uM')

    # Save and close the file so the program can run unattended
    fig.subplots_adjust(hspace=.5)
    plt.savefig(file_name)
    file.close()

    return
