import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib
matplotlib.style.use('axel_style')


def fractured_particles(data_directory):
    data = []
    time = np.array([])
    with open(data_directory + '/fractured_particles.dou', 'r') as file:
        for line in file:
            row = line.strip().split(',')
            row = np.array([float(x) for x in row])  # Create a NumPy array for each row
            data.append(row)
            time= np.append(time, row[-1])
    fracture_array = np.array(data, dtype=object)
    number_of_fractures = [np.shape(x)[0]-1 for x in fracture_array]
    return time, number_of_fractures, fracture_array


if __name__=='__main__':
    # simulation_directory = 'c:/Users/Axel/Documents/DEM/results/fracturing_electrode/SN_1/fracturing_periodic_packing/'
    # simulation_directory = 'c:/Users/Axel/Documents/DEM/results/fracturing_electrode/SN_1/fracturing_electrode_calendering/'
    simulation_directory = 'c:/Users/Axel/Documents/DEM/results/fracturing_electrode/SN_2/fracturing_electrode_calendering/'

    time, number_of_fractures, fracture_array = fractured_particles(simulation_directory)

    figure_fracture, ax_fracture = plt.subplots()
    ax_fracture.plot(time, number_of_fractures)
    ax_fracture.set_xlabel('Time [s]')
    ax_fracture.set_ylabel('Number of fractured particles [-]')


    plt.show()