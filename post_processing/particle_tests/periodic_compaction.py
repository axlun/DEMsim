#
#  Created by Axel on 2025-01-15
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


def box_volume(data_directory):
    periodic_data = np.genfromtxt(data_directory + '/periodic_bc.dou', delimiter=', ')
    surface_area = (periodic_data[:, 2] - periodic_data[:, 1]) * (periodic_data[:, 4] - periodic_data[:,3])
    surface_data = np.genfromtxt(data_directory + '/surface_positions.dou', delimiter=', ')
    box_height = surface_data[:, 5] - surface_data[:, 20]
    box_volume = box_height * surface_area
    return box_volume


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
    ke_data = np.genfromtxt(data_directory + '/kinetic_energy.dou', delimiter=', ')
    time = ke_data[:,-1]
    return time

def kinetic_energy(data_directory):
    ke_data = np.genfromtxt(data_directory + '/kinetic_energy.dou', delimiter=', ')
    ke = ke_data[:, 0]
    return ke

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

def top_surface_pressure(data_directory):
    force_data = np.genfromtxt(data_directory + '/surface_forces.dou', delimiter=', ')
    top_surface_force = force_data[:,1]
    periodic_data = np.genfromtxt(data_directory + '/periodic_bc.dou', delimiter=', ')
    surface_area = (periodic_data[:, 2] - periodic_data[:, 1]) * (periodic_data[:, 4] - periodic_data[:,3])
    top_surface_pressure = top_surface_force / surface_area
    return top_surface_pressure

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
    volume = np.asarray(volume)
    return volume


def fractured_particles(data_directory):
    data = []
    with open(data_directory + '/fractured_particles.dou', 'r') as file:
        for line in file:
            row = line.strip().split(',')
            row = np.array([float(x) for x in row])  # Create a NumPy array for each row
            data.append(row)
    array = np.array(data, dtype=object)
    return array


def relative_density(data_directory):
    particles_vol = particles_volume(data_directory)
    box_v = box_volume(data_directory)
    rel_density = particles_vol / box_v
    return rel_density

def data_grabber(data_directory):
    data_points = np.genfromtxt(data_directory, delimiter=',	')
    ydata = data_points[:, 1]
    xdata = data_points[:,0]
    return xdata, ydata



def main():
    # simulation_directory = "c:/Users/Axel/Documents/DEM/results/particle_tests/fracture_particle_periodic_compaction/SN_0/"
    simulation_directory = "c:/Users/Axel/Documents/DEM/results/particle_tests/fracture_particle_periodic_compaction/SN_1/"
    # simulation_directory = "c:/Users/Axel/Documents/DEM/results/particle_tests/fracture_particle_periodic_compaction/SN_2/"

    box_vol = box_volume(simulation_directory)
    ke = kinetic_energy(simulation_directory)
    particle_vol = particles_volume(simulation_directory)
    pressure = top_surface_pressure(simulation_directory)
    time = Time(simulation_directory)
    rel_density = relative_density(simulation_directory)
    frac_particles = fractured_particles(simulation_directory)
    frac_particles_len = [np.shape(x)[0]-1 for x in frac_particles]


    figure_frac_particles, ax_frac_particles = plt.subplots()
    ax_frac_particles.plot(time, frac_particles_len)
    ax_frac_particles.set_xlabel("Time [s]")
    ax_frac_particles.set_ylabel("Number of fractured particles [-]")

    figure_frac_particles_density, ax_frac_particles_density = plt.subplots()
    ax_frac_particles_density.plot(rel_density, frac_particles_len)
    ax_frac_particles_density.set_xlabel("relative density [-]")
    ax_frac_particles_density.set_ylabel("Number of fractured particles [-]")

    figure_density_time, ax_density_time = plt.subplots()
    ax_density_time.plot(time, rel_density)
    ax_density_time.set_xlabel("Time [s]")
    ax_density_time.set_ylabel("Density [-]")

    figure_volume_time, ax_volume_time = plt.subplots()
    ax_volume_time.plot(time, box_vol)
    ax_volume_time.set_xlabel("Time [s]")
    ax_volume_time.set_ylabel("Volume [m^3]")

    figure_pressure , ax_pressure = plt.subplots()
    ax_pressure.plot(rel_density,pressure, 'r-*',linewidth=3)
    ax_pressure.set_ylabel("Pressure")
    ax_pressure.set_xlabel("Relative density")
    # ax_pressure.legend()

    figure_pressure_time , ax_pressure_time = plt.subplots()
    ax_pressure_time.plot(time ,pressure, 'r-',linewidth=3)
    ax_pressure_time.set_ylabel("Pressure")
    ax_pressure_time.set_xlabel("Time")

    figure_volume, ax_volume = plt.subplots()
    ax_volume.plot(time ,particle_vol, 'r-*',linewidth=3)
    ax_volume.set_ylabel("Particle volume")
    ax_volume.set_xlabel("Time")
    # ax_volume.legend()


    figure_ke, ax_ke = plt.subplots()
    ax_ke.plot(time, ke)
    ax_ke.set_ylabel("Kinetic energy")
    ax_ke.set_xlabel("Time [s]")


    plt.show()

if __name__ == '__main__':
    main()
