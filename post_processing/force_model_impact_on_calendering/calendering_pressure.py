import glob

from math import pi

import numpy as np

import matplotlib.pyplot as plt
import matplotlib

matplotlib.style.use('classic')








if __name__ == '__main__':
    simulation_directory = 'c:/Users/Axel/Documents/DEM/results/electrode_calendaring/bt065N500'


    force_data = np.genfromtxt(simulation_directory + '/surface_forces.dou', delimiter=', ')
    surface_position_data = np.genfromtxt(simulation_directory + '/surface_positions.dou', delimiter=', ')
    with open(simulation_directory + '/surface_positions.dou', 'r') as position_data_file:
        first_line = position_data_file.readlines()[0]
        first_line = first_line.split(', ')
        id_idx = [i for i in range(len(first_line)) if first_line[i].startswith('ID')]
    surface_types = [first_line[idx+1][5:] for idx in id_idx]

    surface_indices = [i for i, surface_type in enumerate(surface_types) if surface_type == 'PointSurface']

    surface1_force = force_data[:, id_idx[0]+1]

    periodic_BC_data = np.genfromtxt(simulation_directory + '/periodic_bc.dou', delimiter=', ')
    periodic_BC_x_min = periodic_BC_data[:,1]
    periodic_BC_x_max = periodic_BC_data[:,2]
    periodic_BC_y_min = periodic_BC_data[:,3]
    periodic_BC_y_max = periodic_BC_data[:,4]
    x_side_leght = periodic_BC_x_max - periodic_BC_x_min
    y_side_leght = periodic_BC_y_max - periodic_BC_y_min
    calendering_pressure = surface1_force/(x_side_leght*y_side_leght)
    calendering_time = periodic_BC_data[:,-1]
    calendering_surface_position = surface_position_data[:,-2]
    print(calendering_pressure)


    plt.plot(calendering_surface_position, calendering_pressure)
    plt.ylabel("Pressure [Pa]")
    plt.xlabel("Calendering surface position [m]")
    plt.show()