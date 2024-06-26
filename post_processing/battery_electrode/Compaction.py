import glob

from math import pi

import numpy as np

import matplotlib.pyplot as plt
import matplotlib

matplotlib.style.use('classic')


def particle_volume():
    p_volume = 0.361901
    return p_volume


def dimensions_box(data_directory):
    with open(data_directory + '/surface_positions.dou', 'r') as position_data_file:
        first_line = position_data_file.readlines()[0]
        first_line = first_line.split(', ')
    id_idx = [i for i in range(len(first_line)) if first_line[i].startswith('ID')]
    surface_types = [first_line[idx+1][5:] for idx in id_idx]
    if surface_types.count('ID=') == 0 and surface_types.count('PointSurface') == 2:
        id_idx.sort(key=lambda x: first_line[x+1])
        wall_data = np.genfromtxt(data_directory + '/surface_positions.dou', delimiter=', ')
        data = np.zeros((wall_data.shape[0], 1))
        data = wall_data[:,  id_idx[0]+32]
        time= wall_data[:,id_idx[0]+42]

        return data
    else:
        raise ValueError("A box could not be defined from the data in " + data_directory +
                         '/surface_positions.dou')


def Time(data_directory):
    with open(data_directory + '/surface_positions.dou', 'r') as position_data_file:
        first_line = position_data_file.readlines()[0]
        first_line = first_line.split(', ')
    id_idx = [i for i in range(len(first_line)) if first_line[i].startswith('ID')]
    surface_types = [first_line[idx+1][5:] for idx in id_idx]
    if surface_types.count('ID=') == 0 and surface_types.count('PointSurface') == 2:
        id_idx.sort(key=lambda x: first_line[x+1])
        wall_data = np.genfromtxt(data_directory + '/surface_positions.dou', delimiter=', ')
        time = np.zeros((wall_data.shape[0], 1))
        time = wall_data[:, id_idx[0]+42]

        return time


def pressures_box(data_directory):
    force_data = np.genfromtxt(data_directory + '/surface_forces.dou', delimiter=', ')
    with open(data_directory + '/surface_positions.dou', 'r') as position_data_file:
        first_line = position_data_file.readlines()[0]
        first_line = first_line.split(', ')
        id_idx = [i for i in range(len(first_line)) if first_line[i].startswith('ID')]
    surface_types = [first_line[idx+1][5:] for idx in id_idx]

    surface_indices = [i for i, surface_type in enumerate(surface_types) if surface_type == 'PointSurface']

    surface1_force = force_data[:, surface_indices[0]+6]
    p = surface1_force/(0.448928*0.448928*4)

    return p


if __name__ == '__main__':
    simulation_directory = 'c:/Users/Axel/Documents/DEM/results/electrode_calendaring/bt065N1500'

    volume_box = (0.448928*2)**2 * dimensions_box(simulation_directory)[:]
    porosity = (1-(particle_volume()*(1+0.07/(0.33+0.07)))/volume_box)
    pressures = pressures_box(simulation_directory)[:]
    time = Time(simulation_directory)[:]
    plt.plot(pressures, porosity*100)
    plt.xlabel("Pressure [Pa]")
    plt.ylabel("Porosity")
    plt.show()
    plt.plot(time, pressures)
    plt.ylabel("Pressure [Pa]")
    plt.xlabel("time")
    plt.show()
# 330-940-1460-1932-2352

