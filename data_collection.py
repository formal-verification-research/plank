from numpy import count_nonzero

def data_collection(model, x_steps, current_time_step, nodes):

    if current_time_step == 0:
        file_data = open("Data", 'w')
    else:
        file_data = open("Data", 'a')

    # Collects the percentage of Extra Cellular Matrix covered by a newly formed capillary
    density = count_nonzero(model[1:, ])
    coverage = density / (x_steps * (x_steps / 2 + .5)) * 100
    file_data.write("Time: " + str(current_time_step * 8 / 3600) + "\n")
    file_data.write("Density:  " + str(density) + "\n")
    file_data.write("Blood vessel coverage (% of area):  " + str(coverage) + "\n")
    file_data.write("Number of nodes: " + str(nodes) + "\n\n")

    file_data.close()

    return
