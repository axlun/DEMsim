#
#  Created by Axel on 2024-02-28
#
import glob
import os

from math import pi

import numpy as np

import matplotlib.pyplot as plt
import matplotlib
import re

import pandas as pd

matplotlib.style.use('axel_style')


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
    particle_files = os.listdir(data_directory + '/particles/')
    time = []
    particle_time_and_file_dict = {}
    for i in particle_files:
        time_stamp = re.split(r'\Aparticles_', i)
        time_stamp = re.split(r'.dou\Z', time_stamp[1])
        time.append(float(time_stamp[0]))
        particle_time_and_file_dict[time_stamp[0]] = i
    time.sort()
    volume = []
    for i in time:
        key = str(i)
        if i.is_integer(): key = str(int(i))

        particle_data = pd.read_csv(
            data_directory + 'particles/' + particle_time_and_file_dict[key]).to_numpy()
        particle_radius = particle_data[:, 7]

        temp_volume = 0
        for i in particle_radius:
            temp_volume += pi * i ** 3 * 4 / 3
        volume.append(temp_volume)
    return volume


def relative_density(data_directory):
    particles_vol = particles_volume(data_directory)
    box_volume = (dimensions_box(data_directory) * 2) ** 3
    rel_density = particles_vol / box_volume
    return rel_density

def data_grabber(data_directory):
    data_points = np.genfromtxt(data_directory, delimiter=',	')
    # print(data_points[:, 0])
    #data_points = [data_points[:, 1], data_points[:, 2]]
    ydata = data_points[:, 1]
    xdata = data_points[:,0]
    return xdata, ydata


def main():
    # simulation_directory = 'c:/Users/Axel/Documents/DEM/results/cube_die_compaction/3'
    simulation_directory = "c:/Users/Axel/Documents/DEM/results/particle_tests/swelling/swelling_cube_die_compaction/"
    particle_vol = particles_volume(simulation_directory)
    pressure = pressures_box(simulation_directory)
    time = Time(simulation_directory)
    rel_density = relative_density(simulation_directory)


    figure_pressure , ax_pressure = plt.subplots()
    ax_pressure.plot(rel_density,pressure, 'r-*',linewidth=3)
    ax_pressure.set_ylabel("Pressure")
    ax_pressure.set_xlabel("Relative density")
    # ax_pressure.legend()

    figure_pressure_time , ax_pressure_time = plt.subplots()
    ax_pressure_time.plot(time ,pressure, 'r-*',linewidth=3)
    ax_pressure_time.set_ylabel("Pressure")
    ax_pressure_time.set_xlabel("Time")

    figure_volume, ax_volume = plt.subplots()
    ax_volume.plot(time ,particle_vol, 'r-*',linewidth=3)
    ax_volume.set_ylabel("Particle volume")
    ax_volume.set_xlabel("Time")
    # ax_volume.legend()


    simulation_directory_unload = "c:/Users/Axel/Documents/DEM/results/particle_tests/swelling/" \
                                  "swelling_cube_die_compaction/unload_D=0.8/"
    particle_vol_unload = particles_volume(simulation_directory_unload)
    pressure_unload  = pressures_box(simulation_directory_unload)
    time_unload  = Time(simulation_directory_unload)
    rel_density_unload  = relative_density(simulation_directory_unload)

    ax_pressure.plot(rel_density_unload,pressure_unload, 'b-*',linewidth=3)

    ax_pressure_time.plot(time_unload ,pressure_unload, 'b-*',linewidth=3)

    ax_volume.plot(time_unload ,particle_vol_unload, 'r-*',linewidth=3)


    plt.show()
    # plt.plot(rel_density, pressure, label='DEM')
    # reference_data_dir = 'c:/Users/Axel/Documents/DEM/results/cube_die_compaction/Data_from_Skrinjar/Isostatic_compaction_monolithic_5000_particles'
    # print(data_grabber(reference_data_dir))
    # plt.plot(data_grabber(reference_data_dir)[0], data_grabber(reference_data_dir)[1],label='Skrinjar')
    # plt.legend()
    # plt.ylabel("Pressure")
    # plt.xlabel("Relative density D")
    # plt.xlim(.6, .8)
    # plt.show()



    # plt.plot(time, pressure)
    # plt.ylabel("Pressure [Pa]")
    # plt.xlabel("time")
    # plt.show()

    # plt.plot(time, particle_vol)
    # plt.ylabel("Volume")
    # plt.xlabel("Time")
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
