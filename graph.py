"""
Description:
The graph file takes the different substrate matrices and the model matrix into a graph for better visual
interpretation of the results
"""


# Imports
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits import mplot3d
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from matplotlib import cm
from parameter_vault import V1, F0, L, DP, x_steps
from math import floor


# Function
def graph(y_substrate, vegf, pedf, fib, pro, x_vector, y_vector, model, current_time_step, total_number_time_steps,
          total_time):

    # Create the height data for the 3D graph
    vegf_z = []
    pedf_z = []
    fib_z = []
    pro_z = []
    x_graph = []
    y_graph = []

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

    # Create the x and y axis for the graphs
    for i in range(len(x_vector)):
        x_graph.append(x_vector[i] / x_steps * L)
    for i in range(len(y_vector)):
        y_graph.append(y_vector[i] / y_substrate * L)

    # Create the overall file
    hours = current_time_step / total_number_time_steps * total_time * L * L / DP
    file_name = "Time = " + str(round(hours, 1)) + " Hours" + ".png"
    file = open(file_name, "a+")
    fig = plt.figure(figsize=plt.figaspect(10))

    # Create the EC color map
    new_paired = cm.get_cmap('Paired', 7)
    ax = fig.add_subplot(5, 1, 1)
    ax.imshow(model, cmap=new_paired)
    ax.title.set_text('Angiogenesis')

    # Create the vegf 3D graph
    ax = fig.add_subplot(5, 1, 2, projection='3d')
    ax.plot_trisurf(x_graph, y_graph, vegf_z, cmap='Greens', edgecolor='none')
    ax.title.set_text('VEGF')
    ax.set_xlabel('mm', labelpad=12)
    ax.set_ylabel('mm', labelpad=12)
    ax.set_zlabel('uM', labelpad=12)

    # Create the pedf 3D graph
    ax = fig.add_subplot(5, 1, 3, projection='3d')
    ax.plot_trisurf(x_graph, y_graph, pedf_z, cmap='Blues', edgecolor='none')
    ax.title.set_text('PEDF')
    ax.set_xlabel('mm', labelpad=12)
    ax.set_ylabel('mm', labelpad=12)
    ax.set_zlabel('uM', labelpad=12)

    # Create the fib 3D graph
    ax = fig.add_subplot(5, 1, 4, projection='3d')
    ax.plot_trisurf(x_graph, y_graph, fib_z, cmap='viridis', edgecolor='none')
    ax.set_zlim3d(0, 0.01)
    ax.title.set_text('Fibronectin')
    ax.set_xlabel('mm', labelpad=12)
    ax.set_ylabel('mm', labelpad=12)
    ax.set_zlabel('uM', labelpad=12)

    # Create the protease 3D graph
    ax = fig.add_subplot(5, 1, 5, projection='3d')
    ax.plot_trisurf(x_graph, y_graph, pro_z, cmap='magma', edgecolor='none')
    ax.title.set_text('Protease')
    ax.set_xlabel('mm', labelpad=12)
    ax.set_ylabel('mm', labelpad=12)
    ax.set_zlabel('uM', labelpad=12)

    # Save and close the file so the program can run unattended
    fig.subplots_adjust(hspace=0.7)
    plt.savefig(file_name)
    file.close()

    return
