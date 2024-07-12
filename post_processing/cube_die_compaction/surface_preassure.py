import glob

from math import pi

import numpy as np

import matplotlib.pyplot as plt
import matplotlib

matplotlib.style.use('classic')


def dimensions_box(data_directory):
    with open(data_directory + '/surface_positions.dou', 'r') as position_data_file:
        first_line = position_data_file.readlines()[0]
        first_line = first_line.split(', ')
    id_idx = [i for i in range(len(first_line)) if first_line[i].startswith('ID')]
    surface_types = [first_line[idx + 1][5:] for idx in id_idx]
    if surface_types.count('ID=') == 0 and surface_types.count('PointSurface') == 6:
        id_idx.sort(key=lambda x: first_line[x + 1])
        wall_data = np.genfromtxt(data_directory + '/surface_positions.dou', delimiter=', ')
        data = np.zeros((wall_data.shape[0], 1))
        data = wall_data[:, id_idx[0] + 5]
        return data
    else:
        raise ValueError("A box could not be defined from the data in " + data_directory +
                         '/surface_positions.dou')


def Time(data_directory):
    with open(data_directory + '/surface_positions.dou', 'r') as position_data_file:
        first_line = position_data_file.readlines()[0]
        first_line = first_line.split(', ')
    id_idx = [i for i in range(len(first_line)) if first_line[i].startswith('ID')]
    surface_types = [first_line[idx + 1][5:] for idx in id_idx]
    if surface_types.count('ID=') == 0 and surface_types.count('PointSurface') == 6:
        id_idx.sort(key=lambda x: first_line[x + 1])
        wall_data = np.genfromtxt(data_directory + '/surface_positions.dou', delimiter=', ')
        time = np.zeros((wall_data.shape[0], 1))
        time = wall_data[:, id_idx[0] + 90]
        return time


def pressures_box(data_directory):
    force_data = np.genfromtxt(data_directory + '/surface_forces.dou', delimiter=', ')
    with open(data_directory + '/surface_positions.dou', 'r') as position_data_file:
        first_line = position_data_file.readlines()[0]
        first_line = first_line.split(', ')
        id_idx = [i for i in range(len(first_line)) if first_line[i].startswith('ID')]
    surface_types = [first_line[idx + 1][5:] for idx in id_idx]

    surface_indices = [i for i, surface_type in enumerate(surface_types) if surface_type == 'PointSurface']
    surface1_force = force_data[:, surface_indices[0] + 6]
    box_surface1_area = (dimensions_box(data_directory) * 2) ** 2
    box_surface1_pressure = surface1_force / box_surface1_area

    return box_surface1_pressure


def surface_forces(data_directory):
    force_data = np.genfromtxt(data_directory + '/surface_forces.dou', delimiter=', ')
    with open(data_directory + '/surface_positions.dou', 'r') as position_data_file:
        first_line = position_data_file.readlines()[0]
        first_line = first_line.split(', ')
        id_idx = [i for i in range(len(first_line)) if first_line[i].startswith('ID')]
    surface_types = [first_line[idx + 1][5:] for idx in id_idx]
    surface_indices = [i for i, surface_type in enumerate(surface_types) if surface_type == 'PointSurface']
    surface1_force_res = force_data[:,
                         surface_indices[0] + 4]  # 4th columns after surface_indices is the resultant force
    return surface1_force_res


def particles_volume(data_directory):
    particle_data = np.genfromtxt(data_directory + '/particles/particles_0.dou', delimiter=', ')
    particle_radius = particle_data[:, 7]
    volume = 0
    for i in range(len(particle_radius)):
        volume += pi * particle_radius[i] ** 3 * 4 / 3
    return volume


def relative_density(data_directory):
    particles_vol = particles_volume(data_directory)
    box_volume = (dimensions_box(data_directory) * 2) ** 3
    rel_density = particles_vol / box_volume
    print(rel_density)
    return rel_density

def data_grabber(data_directory):
    data_points = np.genfromtxt(data_directory, delimiter=',	')
    print(data_points[:, 0])
    #data_points = [data_points[:, 1], data_points[:, 2]]
    ydata = data_points[:, 1]
    xdata = data_points[:,0]
    return xdata, ydata


def main():
    # simulation_directory = 'c:/Users/Axel/Documents/DEM/results/cube_die_compaction/3'
    simulation_directory = "C:/Users/Axel/Documents/DEM/results/particle_tests/swelling/swelling_cube_die_compaction/"
    particles_volume(simulation_directory)
    pressure = pressures_box(simulation_directory)
    # time = Time(simulation_directory)
    rel_density = relative_density(simulation_directory)
    plt.plot(rel_density, pressure / 200e6,label='DEM')
    reference_data_dir = 'c:/Users/Axel/Documents/DEM/results/cube_die_compaction/Data_from_Skrinjar/Isostatic_compaction_monolithic_5000_particles'
    print(data_grabber(reference_data_dir))
    plt.plot(data_grabber(reference_data_dir)[0], data_grabber(reference_data_dir)[1],label='Skrinjar')
    plt.legend()
    plt.ylabel("Pressure/Yield strength")
    plt.xlabel("Relative density D")
    plt.xlim(.6, .8)
    plt.show()



    # plt.plot(time, pressure)
    # plt.ylabel("Pressure [Pa]")
    # plt.xlabel("time")
    # plt.show()
    # rint(surface_forces(simulation_directory))
    # volume_box = (0.448928*2)**2 * dimensions_box(simulation_directory)[:]
    # porosity = (1-(particle_volume()*(1+0.07/(0.33+0.07)))/volume_box)
    # pressures = pressures_box(simulation_directory)[:]
    # time = Time(simulation_directory)[:]
    # plt.plot(pressures, porosity*100)
    # plt.xlabel("Pressure [Pa]")
    # plt.ylabel("Porosity")
    # plt.show()
    # plt.plot(time, pressures)
    # plt.ylabel("Pressure [Pa]")
    # plt.xlabel("time")
    # plt.show()


if __name__ == '__main__':
    main()
