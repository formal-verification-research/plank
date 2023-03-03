from numpy import count_nonzero

def data_collection(model, x_steps, current_time_step, nodes, pro, vegf_old, pedf_old, divide_prob):

    if current_time_step == 0:
        file_data = open("Data", 'w')
        file_density = open("Density", 'w')
        file_node = open("Nodes", 'w')
        file_densitynodes = open("DensityVSNodes", 'w')
    else:
        file_data = open("Data", 'a')
        file_density = open("Density", 'a')
        file_node = open("Nodes", 'a')
        file_densitynodes = open("DensityVSNodes", 'a')

    # Collects the percentage of Extra Cellular Matrix covered by a newly formed capillary
    density = count_nonzero(model[1:, ])
    coverage = density / (x_steps * (x_steps / 2 + .5)) * 100
    file_data.write("Time" + "\t" +str(current_time_step * 8 / 3600) + "\n")
    file_data.write("Density" + "\t" + str(density) + "\n")
    file_data.write("Blood vessel coverage (% of area)" + "\t" + str(coverage) + "\n")
    file_data.write("Number of nodes" + "\t" + str(nodes) + "\n\n")

    # Separate collected data into different files for easier anaylsis, First number= time step, Second=Density
    file_density.write(str(current_time_step * 8 / 3600) + "\t" + str(coverage) + "\n")

    # First number = time step, Second = Nodes, Third = VEGF concentration at division, Fourth = Probability of division
    file_node.write(str(current_time_step * 8 / 3600) + "\t" + str(nodes) + "\n")

    file_densitynodes.write(str(current_time_step * 8 / 3600) + "\t" + str(coverage) + "\t" + str(nodes) + "\n")

    file_data.close()

    return
