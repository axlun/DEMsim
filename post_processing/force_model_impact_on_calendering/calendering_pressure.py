from Bertil_functions.Bertil_functions import *

import numpy as np


import matplotlib.pyplot as plt
import matplotlib

matplotlib.style.use('classic')


def local_data_gatherer(simulation_directory):
    periodic_BC_data = np.genfromtxt(simulation_directory + '/periodic_bc.dou', delimiter=', ')
    force_data = np.genfromtxt(simulation_directory + '/surface_forces.dou', delimiter=', ')
    with open(simulation_directory + '/surface_forces.dou', 'r') as force_data_file:
        first_line = force_data_file.readlines()[0]
        first_line = first_line.split(', ')
        id_idx = [i for i in range(len(first_line)) if first_line[i].startswith('ID')]
        surface_force_index = id_idx
    force_data_file.close()
    surface_position_data = np.genfromtxt(simulation_directory + '/surface_positions.dou', delimiter=', ')
    with open(simulation_directory + '/surface_positions.dou', 'r') as position_data_file:
        first_line = position_data_file.readlines()[0]
        first_line = first_line.split(', ')
        id_idx = [i for i in range(len(first_line)) if first_line[i].startswith('ID')]
        surface_position_index = id_idx
    position_data_file.close()
    surface_types = [first_line[idx+1][5:] for idx in id_idx]

    surface_indices = [i for i, surface_type in enumerate(surface_types) if surface_type == 'PointSurface']

    return force_data, surface_force_index, surface_position_index,surface_position_data, periodic_BC_data
def bertil_data_gatherer(simulation_directory):
    ##================SURFACE FORCE DATA===================================================================================
    force_data_list = one_file_reader(simulation_directory + '/surface_forces.dou')
    force_data_list = np.asarray(force_data_list)
    force_data_list = np.char.split(force_data_list, ', ')
    first_line = force_data_list[0]
    id_idx = [i for i in range(len(first_line)) if first_line[i].startswith('ID')]
    surface_force_index = id_idx
    for i in force_data_list[:]:
        if 'force_data' in locals():
            force_data = np.vstack((force_data, np.asarray(i)))
        else:
            force_data = np.asarray(i)
    for i in enumerate(force_data): force_data[i[0], -1] = force_data[i[0], -1][:-1]
    ##================PERIODIC BC DATA==================================================================================
    periodic_BC_data_list = one_file_reader(simulation_directory + '/periodic_bc.dou')
    periodic_BC_data_list = np.asarray(periodic_BC_data_list)
    periodic_BC_data_list = np.char.split(periodic_BC_data_list, ', ')
    for i in periodic_BC_data_list[:]:
        if 'periodic_BC_data' in locals():
            periodic_BC_data = np.vstack((periodic_BC_data, np.asarray(i)))
        else:
            periodic_BC_data = np.asarray(i)
    for i in enumerate(periodic_BC_data): periodic_BC_data[i[0], -1] = periodic_BC_data[i[0], -1][:-1]
    periodic_BC_data = periodic_BC_data.astype(float)

    # ================SURFACE POSITION DATA==========================================================================
    surface_position_data_list = one_file_reader(simulation_directory + '/surface_positions.dou')

    first_line = surface_position_data_list[0].split(', ')
    id_idx = [i for i in range(len(first_line)) if first_line[i].startswith('ID')]
    surface_position_index = id_idx
    surface_types = [first_line[idx + 1][5:] for idx in id_idx]
    surface_indices = [i for i, surface_type in enumerate(surface_types) if surface_type == 'PointSurface']

    surface_position_data_list = np.asarray(surface_position_data_list)
    surface_position_data_list = np.char.split(surface_position_data_list, ', ')
    for i in surface_position_data_list[:]:
        if 'surface_position_data' in locals():
            surface_position_data = np.vstack((surface_position_data, np.asarray(i)))
        else:
            surface_position_data = np.asarray(i)
    for i in enumerate(surface_position_data):
        surface_position_data[i[0], -1] = surface_position_data[i[0], -1][:-1]


    return force_data, surface_force_index, surface_position_index,surface_position_data, periodic_BC_data


if __name__ == '__main__':
    #simulation_directory = 'c:/Users/Axel/Documents/DEM/results/electrode_calendaring/SN00'
#    simulation_directory = '/scratch/users/axlun/DEMsim/results/electrode_calendaring/bt065N3000bt_coeff_1_mu_0_mu_wall_0_fraction_binder_contacts_06'
    simulation_directory = '/scratch/users/axlun/DEMsim/results/electrode_calendaring/SN00_1000p_plastic_binder'
#    simulation_directory = 'c:/Users/Axel/Documents/DEM/results/electrode_calendaring/SN0_1000particlesWithBinder_new_packingmethod_2'

    if simulation_directory.startswith("/scratch"):force_data, surface_force_index, surface_position_index, surface_position_data, periodic_BC_data = bertil_data_gatherer(simulation_directory)
    elif simulation_directory.startswith("c:"):force_data, surface_force_index, surface_position_index, surface_position_data, periodic_BC_data = local_data_gatherer(simulation_directory)
    else: print("Error with simulation directory")

    print(surface_force_index[1])
    calendering_surface_force = force_data[:, surface_force_index[1]+1].astype(float)
    periodic_BC_x_min = periodic_BC_data[:,1].astype(float)
    periodic_BC_x_max = periodic_BC_data[:,2].astype(float)
    periodic_BC_y_min = periodic_BC_data[:,3].astype(float)
    periodic_BC_y_max = periodic_BC_data[:,4].astype(float)
    x_side_leght = periodic_BC_x_max - periodic_BC_x_min
    y_side_leght = periodic_BC_y_max - periodic_BC_y_min
    calendering_surface_pressure = calendering_surface_force/(x_side_leght*y_side_leght)
    calendering_time = periodic_BC_data[:,0]
    calendering_surface_position = surface_position_data[:,surface_position_index[1]+14].astype(float)

    fig, ax1 = plt.subplots()
    ax1.set_ylabel("Calendering surface pressure [MPa]")
    ax1.set_xlabel("time [s]")
    ax1.plot(calendering_time, calendering_surface_pressure*1e-6)
    ax2 = ax1.twinx()
#    ax2.set_xlabel("Calendering surface position [m]")
    ax2.set_ylabel("Calendering surface position [m]")
    ax2.plot(calendering_time, calendering_surface_position)
    # ax1.title('calendering surface pressure')
    fig.tight_layout()
    plt.show()


    plt.plot(calendering_surface_position, calendering_surface_pressure*1e-6)
    plt.ylabel("Calendering surface pressure [MPa]")
    plt.xlabel("Calendering surface position [m]")
    plt.title("calendering surface pressure")
    plt.show()