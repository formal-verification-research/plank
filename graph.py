# Description
# graph turns the different substrate matrices and the workspace matrix into a graph for better visual 
# interpretation of the results


# Imports
from matplotlib import pyplot
from mpl_toolkits import mplot3d
from matplotlib import cm


# Function
def graph(y_substrate, x_steps, vegf, pedf, fib, pro, x_vector, y_vector, workspace, current_time_step):

    VEGFzVector = []  # vegf height data for the 3D graph
    PEDFzVector = []  # pedf height data for the 3D graph
    FIBzVector = []  # fib height data for the 3D graph
    PROzVector = []  # pro height data for the 3D graph

    # Add the appropriate info for odd and even lines to the 3D vector
    for y in range(y_substrate):
        if y % 2 == 0:
            for x in range(x_steps-1):
                VEGFzVector.append(vegf[y][x])
                PEDFzVector.append(pedf[y][x])
                FIBzVector.append(fib[y][x])
                PROzVector.append(pro[y][x])
        else:
            for x in range(x_steps):
                VEGFzVector.append(vegf[y][x])
                PEDFzVector.append(pedf[y][x])
                FIBzVector.append(fib[y][x])
                PROzVector.append(pro[y][x])

    # Create the overall file
    file_name = "Time = " + str(current_time_step) + ".pdf"
    file = open(file_name, "a+")
    fig = pyplot.figure(figsize=pyplot.figaspect(4.0))

    # Create the EC color map
    ax = fig.add_subplot(5, 1, 1)
    ax.imshow(workspace, cmap='Purples')
    ax.title.set_text('ECM')

    # Create the vegf 3D graph
    ax = fig.add_subplot(5, 1, 2, projection='3d')
    ax.plot_trisurf(x_vector, y_vector, VEGFzVector, cmap='hot', edgecolor='none')
    ax.set_zlim3d(0, 10)
    ax.title.set_text('VEGF')

    # Create the pedf 3D graph
    ax = fig.add_subplot(5, 1, 3, projection='3d')
    ax.plot_trisurf(x_vector, y_vector, PEDFzVector, cmap='cool', edgecolor='none')
    ax.set_zlim3d(0, 10)
    ax.title.set_text('PEDF')

    # Create the fib 3D graph
    ax = fig.add_subplot(5, 1, 4, projection='3d')
    ax.plot_trisurf(x_vector, y_vector, FIBzVector, cmap='plasma', edgecolor='none')
    ax.set_zlim3d(0, 1)
    ax.title.set_text('Fibronectin')

    # Create the protease 3D graph
    ax = fig.add_subplot(5, 1, 5, projection='3d')
    ax.plot_trisurf(x_vector, y_vector, PROzVector, cmap='viridis', edgecolor='none')
    ax.set_zlim3d(0, 10)
    ax.title.set_text('Protease')

    # Save and close the file so the program can run unattended
    fig.subplots_adjust(hspace=.5)
    pyplot.savefig(file_name)
    file.close()

    return
