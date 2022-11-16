from Bertil_functions.Bertil_functions import *

import numpy as np


import matplotlib.pyplot as plt
import matplotlib

matplotlib.style.use('classic')


def local_data_gatherer(simulation_directory):
    periodic_BC_data = np.genfromtxt(simulation_directory + '/periodic_bc.dou', delimiter=', ')
    force_data = np.genfromtxt(simulation_directory + '/surface_forces.dou', delimiter=', ')
    force_fabric_tensor_data = np.genfromtxt(simulation_directory + '/force_fabric_tensor.dou', delimiter=',')
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



    return force_data, surface_force_index, surface_position_index,surface_position_data, periodic_BC_data, force_fabric_tensor_data


def calendering_break_index_func(calendering_surface_position):
    calendering_initiate_index = 0
    calendering_break_index = -1
    h_al = 1.11 # Active layer height
    flag1 = True
    flag2 = True
    val_hist = -1
    val_hist2 = -1
#    calendering_break_index = 0
    for count,val in enumerate(calendering_surface_position):
       # print(count,val)
        if val < h_al*1.2 and flag1:
            calendering_initiate_index = count
            flag1 = False
        if (val == val_hist2) and (flag1 == False) and flag2:
            flag2 = False
            calendering_break_index = count
        val_hist2= val_hist
        val_hist = val
    return calendering_initiate_index,calendering_break_index

def calendering_surface_force_func(force_data,surface_force_index):
    calendering_surface_force = force_data[:, surface_force_index[1]+1].astype(float)
    return calendering_surface_force

def calendering_surface_pressure_func(periodic_BC_data,force_data,surface_force_index,surface_position_index):
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
    return calendering_time,calendering_surface_pressure,calendering_surface_position

if __name__ == '__main__':

    simulation_directory = 'c:/Users/Axel/Documents/DEM/results/natural_packing_hertz/cm_particles_N_200_mass_sclaing_1e1_gravity_1e1_time_step_1e0_run_2/'

    force_data = np.genfromtxt(simulation_directory + '/surface_forces.dou', delimiter=', ')
    force_fabric_tensor_data = np.genfromtxt(simulation_directory + '/force_fabric_tensor.dou', delimiter=',')
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

    # ##================SURFACE FORCE DATA===================================================================================
    # force_data_list = one_file_reader(simulation_directory + '/surface_forces.dou')
    # force_data_list = np.asarray(force_data_list)
    # force_data_list = np.char.split(force_data_list, ', ')
    # first_line = force_data_list[0]
    # id_idx_surf_force = [i for i in range(len(first_line)) if first_line[i].startswith('ID')]
    # surface_force_index = id_idx_surf_force
    # for i in force_data_list[:]:
    #     if 'force_data' in locals():
    #         force_data = np.vstack((force_data, np.asarray(i)))
    #     else:
    #         force_data = np.asarray(i)
    # for i in enumerate(force_data): force_data[i[0], -1] = force_data[i[0], -1][:-1]
    #
    # # ================SURFACE POSITION DATA==========================================================================
    # surface_position_data_list = one_file_reader(simulation_directory + '/surface_positions.dou')
    #
    # first_line = surface_position_data_list[0].split(', ')
    # id_idx_suf_pos = [i for i in range(len(first_line)) if first_line[i].startswith('ID')]
    # surface_position_index = id_idx_suf_pos
    # surface_types = [first_line[idx + 1][5:] for idx in id_idx_suf_pos]
    # surface_indices = [i for i, surface_type in enumerate(surface_types) if surface_type == 'PointSurface']
    #
    # surface_position_data_list = np.asarray(surface_position_data_list)
    # surface_position_data_list = np.char.split(surface_position_data_list, ', ')
    # for i in surface_position_data_list[:]:
    #     if 'surface_position_data' in locals():
    #         surface_position_data = np.vstack((surface_position_data, np.asarray(i)))
    #     else:
    #         surface_position_data = np.asarray(i)
    # for i in enumerate(surface_position_data):
    #     surface_position_data[i[0], -1] = surface_position_data[i[0], -1][:-1]
    # ##================FORCE FABRIC TENSOR DATA==================================================================================
    # force_fabric_tensor_data_list = one_file_reader(simulation_directory + '/force_fabric_tensor.dou')
    # force_fabric_tensor_data_list = np.asarray(force_fabric_tensor_data_list)
    # force_fabric_tensor_data_list = np.char.split(force_fabric_tensor_data_list, ', ')
    # for i in force_fabric_tensor_data_list[:]:
    #     if 'force_fabric_tensor_data' in locals():
    #         force_fabric_tensor_data = np.vstack((force_fabric_tensor_data, np.asarray(i)))
    #     else:
    #         force_fabric_tensor_data = np.asarray(i)
    # for i in enumerate(force_fabric_tensor_data): force_fabric_tensor_data[i[0], -1] = force_fabric_tensor_data[i[0], -1][:-1]
    # force_fabric_tensor_data = force_fabric_tensor_data.astype(float)

#    print(surface_force_index[1])
    calendering_surface_force = force_data[:, surface_force_index[1]+1].astype(float)
    bottom_surface_force = force_data[:, surface_force_index[0]+1].astype(float)
    x_side_length = abs(2*surface_position_data[:,surface_position_index[2]+3].astype(float))
    y_side_length = x_side_length

    calendering_surface_pressure = calendering_surface_force/(x_side_length*y_side_length)
    bottom_surface_pressure = bottom_surface_force/(x_side_length*y_side_length)
    calendering_time = surface_position_data[:,-1].astype(float)
    calendering_surface_position = surface_position_data[:,surface_position_index[1]+14].astype(float)
#=================================STRESS IN Z-DIRECTION===============================================================

    vol = x_side_length * y_side_length * calendering_surface_position
    szz = -force_fabric_tensor_data[:, 9] / vol
    # print(vol)
    # print(szz)

    # plt.ion


    fig_calendering_surface_bottom_surface,ax_calendering_surface_bottom_surface = plt.subplots()
    lns_calendering_surface_pressure_2 = ax_calendering_surface_bottom_surface.plot(calendering_time, calendering_surface_pressure * 1e-6,'g', linewidth=2, label=r'Surface pressure')
    lns_bottom_surface_pressure = ax_calendering_surface_bottom_surface.plot(calendering_time, bottom_surface_pressure * 1e-6,'b', linewidth=2 ,label=r'Bottom pressure')
    lns_macroscopic_stress = ax_calendering_surface_bottom_surface.plot(calendering_time, -szz * 1e-6, 'y', linewidth=2, label=r'$\sigma_{zz}$')

    ax_calendering_surface_position_2 = ax_calendering_surface_bottom_surface.twinx()
    ax_calendering_surface_position_2.set_ylabel("Calendering surface position [µm]")
#    ax_calendering_surface_position_2.set_ylim([0,200])

    lns_calendering_surface_position_2 = ax_calendering_surface_position_2.plot(calendering_time, 1e2*calendering_surface_position, 'r', linewidth=2 , label='Surface position')
    ax_calendering_surface_bottom_surface.set_ylabel("Stress in z [MPa]")
    ax_calendering_surface_bottom_surface.set_xlabel("Time [s]")
    ax_calendering_surface_bottom_surface.set_title("Pressure on calendering surface and bottom surface")
    lns_sum = lns_calendering_surface_pressure_2 + lns_bottom_surface_pressure + lns_macroscopic_stress + lns_calendering_surface_position_2
    labs4 = [l.get_label() for l in lns_sum]
    ax_calendering_surface_bottom_surface.legend(lns_sum, labs4, loc="best")


    plt.show()