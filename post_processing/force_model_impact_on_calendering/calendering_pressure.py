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
#    simulation_directory = 'c:/Users/Axel/Documents/DEM/results/electrode_calendaring/SN00'
#    simulation_directory = '/scratch/users/axlun/DEMsim/results/electrode_calendaring/SN00_3000p_plastic_binder_btr_075_hal_10_bstiffnesscoeff_15'
#    simulation_directory = 'c:/Users/Axel/Documents/DEM/results/electrode_calendaring_hertz/SN00_500particles_hertz_particle_contact'
#    simulation_directory = 'c:/Users/Axel/Documents/DEM/results/electrode_calendaring/SN0_1000particlesWithBinder_new_packingmethod_2'
#    simulation_directory = 'c:/Users/Axel/Desktop/SN00_1500p_plastic_binder_btr_1_hal_08_run_2'
#    simulation_directory = '/scratch/users/axlun/DEMsim/results/electrode_calendaring_restart/SN00_3000p_plastic_binder_btr_065_hal_105'
    simulation_directory = '/scratch/users/axlun/DEMsim/results/electrode_calendaring/SN_hertz_3000p_btr_065_new_Ft_b'
#    simulation_directory = '/scratch/users/axlun/DEMsim/results/electrode_calendaring_restart/SN_hertz_5000p_btr_00_hal_295'

#   simulation_directory='c:/Users/Axel/Desktop/electrode_calendaring_restart/SN00_1500p_plastic_binder_btr_08_hal_10'

#    simulation_directory = '/scratch/users/axlun/DEMsim/results/electrode_calendaring/SN00_1500p_plastic_binder_btr_08'#_hal_08_run_2'

    if simulation_directory.startswith("/scratch"): force_data, surface_force_index, surface_position_index, surface_position_data, periodic_BC_data,force_fabric_tensor_data = bertil_data_gatherer(simulation_directory)
    elif simulation_directory.startswith("c:"): force_data, surface_force_index, surface_position_index, surface_position_data, periodic_BC_data,force_fabric_tensor_data = local_data_gatherer(simulation_directory)
    else: print("Error with simulation directory")

#    print(surface_force_index[1])
    calendering_surface_force = force_data[:, surface_force_index[1]+1].astype(float)
    periodic_BC_x_min = periodic_BC_data[:,1].astype(float)
    periodic_BC_x_max = periodic_BC_data[:,2].astype(float)
    periodic_BC_y_min = periodic_BC_data[:,3].astype(float)
    periodic_BC_y_max = periodic_BC_data[:,4].astype(float)
    x_side_length = periodic_BC_x_max - periodic_BC_x_min
    y_side_length = periodic_BC_y_max - periodic_BC_y_min
    calendering_surface_pressure = calendering_surface_force/(x_side_length*y_side_length)
    calendering_time = periodic_BC_data[:,0]
    calendering_surface_position = surface_position_data[:,surface_position_index[1]+14].astype(float)

#=================================STRESS IN Z-DIRECTION===============================================================

    vol = x_side_length * y_side_length * calendering_surface_position
    szz = -force_fabric_tensor_data[:, 9] / vol
    # print(vol)
    # print(szz)

    # plt.ion
    fig_calendering_surface_pressure, ax_calendering_surface_pressure = plt.subplots()
    ax_calendering_surface_pressure.set_ylabel("Calendering surface pressure [MPa]")
    ax_calendering_surface_pressure.set_xlabel("time [s]")
    lns_calendering_surface_pressure = ax_calendering_surface_pressure.plot(calendering_time, calendering_surface_pressure * 1e-6, 'r', label='Pressure')
    ax_calendering_surface_position = ax_calendering_surface_pressure.twinx()
    ax_calendering_surface_position.set_ylabel("Calendering surface position [m]")
    lns_calendering_surface_position = ax_calendering_surface_position.plot(calendering_time, calendering_surface_position,'b',label='Position')
    ax_calendering_surface_pressure.set_title('calendering surface pressure')

    lns = lns_calendering_surface_pressure + lns_calendering_surface_position
    labs = [l.get_label() for l in lns]
    ax_calendering_surface_pressure.legend(lns, labs, loc='best')
    fig_calendering_surface_pressure.tight_layout()

    calendering_initiate_index, calendering_break_index = calendering_break_index_func(calendering_surface_position)
    fig_calendering_process,ax_calendering_process = plt.subplots()
    ax_calendering_process.plot(1e2*calendering_surface_position[calendering_initiate_index:(calendering_break_index-1)], calendering_surface_pressure[calendering_initiate_index:(calendering_break_index-1)]*1e-6,linewidth=3)
    ax_calendering_process.set_ylabel("Calendering surface pressure [MPa]")
    ax_calendering_process.set_xlabel("Calendering surface position [µm]")
    ax_calendering_process.set_title('calendering surface pressure')
    fig_calendering_process.tight_layout


    fig_full_calendering_sim, ax_full_calendering_sim = plt.subplots()
    ax_full_calendering_sim.plot(calendering_surface_position[:], calendering_surface_pressure[:] * 1e-6)
    ax_full_calendering_sim.set_ylabel("Calendering surface pressure [MPa]")
    ax_full_calendering_sim.set_xlabel("Calendering surface position [m]")
    ax_full_calendering_sim.set_title('calendering surface pressure')

    fig_sigma_zz_calendering_surface_position,ax_sigma_zz_calendering_surface_position = plt.subplots()
    ax_sigma_zz_calendering_surface_position.plot(calendering_surface_position[:],szz)
    ax_sigma_zz_calendering_surface_position.set_ylabel("Stress in z [MPa]")
    ax_sigma_zz_calendering_surface_position.set_xlabel("Calendering surface position [m]")
    ax_sigma_zz_calendering_surface_position.set_title("Stress in Z to calendering surface position")

    fig_sigma_zz_time,ax_sigma_zz_time = plt.subplots()
    lns_sigma_zz_time = ax_sigma_zz_time.plot(calendering_time, szz*1e-6, 'r',label=r'$\sigma_{zz}$')
    lsn_cal_surf_pre_time = ax_sigma_zz_time.plot(calendering_time, calendering_surface_pressure * 1e-6,'g',label=r'Surface pressure')
    ax_sigma_zz_time.set_ylabel("Stress in z [MPa]")
    ax_sigma_zz_time.set_xlabel("time [s]")
    ax_calendering_surface_position_time = ax_sigma_zz_time.twinx()
    ax_calendering_surface_position_time.set_ylabel("Calendering surface position [m]")
    lns_calendering_surface_position_time = ax_calendering_surface_position_time.plot(calendering_time, calendering_surface_position, 'b', label='Surface position')
    lns2 = lns_sigma_zz_time+ lsn_cal_surf_pre_time + lns_calendering_surface_position_time
    labs2 = [l.get_label() for l in lns2]
    ax_sigma_zz_time.legend(lns2, labs2, loc="best")

    fig_calendering_surface_pressure_sigma_zz,ax_sigma_zz_time2 = plt.subplots()
    lns5 = ax_sigma_zz_time2.plot(calendering_time, szz*1e-6, 'r',linewidth=3,label=r'$\sigma_{zz}$')
    ax_sigma_zz_time2.set_ylabel("Stress [MPa]")
    ax_sigma_zz_time2.set_xlabel("time [s]")
    lns6 = ax_sigma_zz_time2.plot(calendering_time, calendering_surface_pressure[:] * 1e-6, 'b',linewidth=3, label='Surface pressure')
  #  ax9.set_title('calendering surface pressure')
    lns3 = lns5 + lns6
    labs3 = [l.get_label() for l in lns3]
    ax_sigma_zz_time2.legend(lns3, labs3, loc="best")

    plt.show()