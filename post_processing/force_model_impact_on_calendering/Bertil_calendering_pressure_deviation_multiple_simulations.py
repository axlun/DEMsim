from Bertil_functions.Bertil_functions import *

import numpy as np
import matplotlib
from scipy import interpolate
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from os.path import exists
import shutil
import os


def calendering_break_index_func(calendering_surface_position):
    calendering_initiate_index = 0
    calendering_break_index = -1
    h_al = 1.11  # Active layer height
    flag1 = True
    flag2 = True
    val_hist = -1
    val_hist2 = -1
    val_hist3 = -1
    #    calendering_break_index = 0
    for count, val in enumerate(calendering_surface_position):
        if val < h_al * 1.2 and flag1:
            calendering_initiate_index = count
            flag1 = False
        if (val == val_hist3) and (flag1 == False) and flag2:
            flag2 = False
            calendering_break_index = count
        val_hist3 = val_hist2
        val_hist2 = val_hist
        val_hist = val
    return calendering_initiate_index, calendering_break_index

def calendering_plot_processing(simulation_directory):
    # =======================================================================================================================
    force_data_SN_run_1, surface_force_index_SN_run_1, surface_position_index_SN_run_1, surface_position_data_SN_run_1, periodic_BC_data_SN_run_1, \
    force_fabric_tensor_data_SN_run_1, kinetic_energy_data = bertil_data_gatherer(
        simulation_directory)

    calendering_surface_force_SN_run_1 = force_data_SN_run_1[:, surface_force_index_SN_run_1[1] + 1].astype(float)
    bottom_surface_force_SN_run_1 = force_data_SN_run_1[:, surface_force_index_SN_run_1[0] + 1].astype(float)
    x_side_length_SN_run_1 = (
                (surface_position_data_SN_run_1[:, surface_position_index_SN_run_1[0] + 9]).astype(float) - (
            surface_position_data_SN_run_1[:, surface_position_index_SN_run_1[0] + 3]).astype(float))
    y_side_length_SN_run_1 = (
                (surface_position_data_SN_run_1[:, surface_position_index_SN_run_1[0] + 10]).astype(float) - (
            surface_position_data_SN_run_1[:, surface_position_index_SN_run_1[0] + 4]).astype(float))

    calendering_surface_pressure = calendering_surface_force_SN_run_1 / (
                x_side_length_SN_run_1 * y_side_length_SN_run_1)
    bottom_surface_pressure = bottom_surface_force_SN_run_1 / (x_side_length_SN_run_1 * y_side_length_SN_run_1)
    time = surface_position_data_SN_run_1[:, -1].astype(float)
    calendering_surface_position = surface_position_data_SN_run_1[:,
                                            surface_position_index_SN_run_1[1] + 14].astype(float)

    # =================================STRESS IN Z-DIRECTION============================================================
    vol_SN_run_1 = x_side_length_SN_run_1 * y_side_length_SN_run_1 * calendering_surface_position
    sig_x = -force_fabric_tensor_data_SN_run_1[:, 1] / vol_SN_run_1
    sig_y = -force_fabric_tensor_data_SN_run_1[:, 5] / vol_SN_run_1
    sig_z = -force_fabric_tensor_data_SN_run_1[:, 9] / vol_SN_run_1
    tau_xy = -force_fabric_tensor_data_SN_run_1[:, 2] / vol_SN_run_1
    tau_xz = -force_fabric_tensor_data_SN_run_1[:, 3] / vol_SN_run_1

    tau_yx = -force_fabric_tensor_data_SN_run_1[:, 4] / vol_SN_run_1
    tau_yz = -force_fabric_tensor_data_SN_run_1[:, 6] / vol_SN_run_1

    tau_zx = -force_fabric_tensor_data_SN_run_1[:, 7] / vol_SN_run_1
    tau_zy = -force_fabric_tensor_data_SN_run_1[:, 8] / vol_SN_run_1

    # =======================================================================================================================


    return time,calendering_surface_pressure,bottom_surface_pressure, calendering_surface_position, sig_x, sig_y,sig_z, tau_xy,tau_xz,tau_yx,tau_yz, tau_zx,tau_zy, kinetic_energy_data


def max_displacement_index_func(surface_position_vec):
    i = 1

    while surface_position_vec[i-1] >= surface_position_vec[i]:
        i += 1

    index_of_max_displacement = i-1
    return index_of_max_displacement


def load_unload_height_func(calendering_surface_pressure,calendering_surface_position):
    load_limit = 0.5*1E6

    load_limit_index = 0
    unload_limit_index = 0

    while calendering_surface_pressure[load_limit_index] < load_limit:
        load_limit_index += 1
        # print(calendering_surface_pressure[load_limit_index])
        unload_limit_index = load_limit_index + 1
    while calendering_surface_pressure[unload_limit_index] > load_limit:
        unload_limit_index += 1

    loading_height = calendering_surface_position[load_limit_index]
    unload_height  = calendering_surface_position[unload_limit_index]
    min_height = min(calendering_surface_position)
    specific_elastic_springback = (unload_height-min_height)/(loading_height-min_height)
    print(loading_height)

    print(min_height)

    print(unload_height)

    return loading_height, unload_height, specific_elastic_springback


if __name__ == '__main__':

    # ==NATUAL PACKING =================================================================================================
    # simulation_directory = '/scratch/users/axlun/DEMsim/results/electrode_natural_packing_hertz/SN_run_1'

    # ==CALENDERING=====================================================================================================
    simulation_directory_SN_run_1 = '/scratch/users/axlun/DEMsim/results/final_runs/SN_run_1/electrode_calendering_hertz'
    simulation_directory_SN_run_1_spread_1 = '/scratch/users/axlun/DEMsim/results/final_runs/SN_run_1_spread_1/electrode_calendering_hertz'
    simulation_directory_SN_run_1_spread_2 = '/scratch/users/axlun/DEMsim/results/final_runs/SN_run_1_spread_2/electrode_calendering_hertz'
    simulation_directory_SN_run_1_spread_3 = '/scratch/users/axlun/DEMsim/results/final_runs/SN_run_1_spread_3/electrode_calendering_hertz'
    simulation_directory_SN_run_1_spread_4 = '/scratch/users/axlun/DEMsim/results/final_runs/SN_run_1_spread_4/electrode_calendering_hertz'
    simulation_directory_SN_run_1_El_Pl = '/scratch/users/axlun/DEMsim/results/final_runs/SN_run_1_El_Pl/electrode_calendering_hertz'
    simulation_directory_SN_run_1_br_175 = '/scratch/users/axlun/DEMsim/results/final_runs/SN_run_1_br_175/electrode_calendering_hertz'
    simulation_directory_SN_run_1_br_10 = '/scratch/users/axlun/DEMsim/results/final_runs/SN_run_1_br_10/electrode_calendering_hertz'
    simulation_directory_SN_run_1_br_05 = '/scratch/users/axlun/DEMsim/results/final_runs/SN_run_1_br_05/electrode_calendering_hertz'
    simulation_directory_SN_run_1_2E = '/scratch/users/axlun/DEMsim/results/final_runs/SN_run_1_2E/electrode_calendering_hertz'

    # ==MECHANICAL LOADING==============================================================================================
    # simulation_directory = '/scratch/users/axlun/DEMsim/results/electrode_mechanical_loading_hertz/SN_hertz_5000p_btr_8_brr_08_dt_1e1_MS_1e2_SR_2e-3_compression'

    # ==RESTING=========================================================================================================
    # simulation_directory = '/scratch/users/axlun/DEMsim/results/electrode_resting_hertz/SN_hertz_5000p_btr_8_brr_08_dt_5e1_MS_1e4_RT_10'

    # ==PLOT PARAMETERS=================================================================================================
    fig_dir = 'C:/temp/figures/Bertil_calendering_pressure_statistical_deviation/'
    try:
        shutil.rmtree(fig_dir)
        os.mkdir(fig_dir)
    except FileNotFoundError:
        print('Directory does not exist')
        os.mkdir(fig_dir)
    except:
        print('Directory could not be removed')
        quit()
    # ==PLOT PARAMETERS=================================================================================================
    plt.style.use('axel_style')

    # ==PLOT DATA=================================================================================================
    calendering_time_SN_run_1, calendering_surface_pressure_SN_run_1, bottom_surface_pressure_SN_run_1, calendering_surface_position_SN_run_1, sxx_SN_run_1, syy_SN_run_1, szz_SN_run_1, tau_xy_SN_run_1, tau_xz_SN_run_1, tau_yx_SN_run_1, tau_yz_SN_run_1, tau_zx_SN_run_1, tau_zy_SN_run_1,kinetic_energy_data_SN_run_1 = calendering_plot_processing(simulation_directory_SN_run_1)

    calendering_time_SN_run_1_spread_1, calendering_surface_pressure_SN_run_1_spread_1, bottom_surface_pressure_SN_run_1_spread_1, calendering_surface_position_SN_run_1_spread_1, sxx_SN_run_1_spread_1, syy_SN_run_1_spread_1, szz_SN_run_1_spread_1, tau_xy_SN_run_1_spread_1, tau_xz_SN_run_1_spread_1, tau_yx_SN_run_1_spread_1, tau_yz_SN_run_1_spread_1, tau_zx_SN_run_1_spread_1, tau_zy_SN_run_1_spread_1, kinetic_energy_data_SN_run_1_spread_1 = calendering_plot_processing(
        simulation_directory_SN_run_1_spread_1)

    calendering_time_SN_run_1_spread_2, calendering_surface_pressure_SN_run_1_spread_2, bottom_surface_pressure_SN_run_1_spread_2, calendering_surface_position_SN_run_1_spread_2, sxx_SN_run_1_spread_2, syy_SN_run_1_spread_2, szz_SN_run_1_spread_2, tau_xy_SN_run_1_spread_2, tau_xz_SN_run_1_spread_2, tau_yx_SN_run_1_spread_2, tau_yz_SN_run_1_spread_2, tau_zx_SN_run_1_spread_2, tau_zy_SN_run_1_spread_2, kinetic_energy_data_SN_run_1_spread_2 = calendering_plot_processing(
        simulation_directory_SN_run_1_spread_2)

    calendering_time_SN_run_1_spread_3, calendering_surface_pressure_SN_run_1_spread_3, bottom_surface_pressure_SN_run_1_spread_3, calendering_surface_position_SN_run_1_spread_3, sxx_SN_run_1_spread_3, syy_SN_run_1_spread_3, szz_SN_run_1_spread_3, tau_xy_SN_run_1_spread_3, tau_xz_SN_run_1_spread_3, tau_yx_SN_run_1_spread_3, tau_yz_SN_run_1_spread_3, tau_zx_SN_run_1_spread_3, tau_zy_SN_run_1_spread_3, kinetic_energy_data_SN_run_1_spread_3 = calendering_plot_processing(
        simulation_directory_SN_run_1_spread_3)

    calendering_time_SN_run_1_spread_4, calendering_surface_pressure_SN_run_1_spread_4, bottom_surface_pressure_SN_run_1_spread_4, calendering_surface_position_SN_run_1_spread_4, sxx_SN_run_1_spread_4, syy_SN_run_1_spread_4, szz_SN_run_1_spread_4, tau_xy_SN_run_1_spread_4, tau_xz_SN_run_1_spread_4, tau_yx_SN_run_1_spread_4, tau_yz_SN_run_1_spread_4, tau_zx_SN_run_1_spread_4, tau_zy_SN_run_1_spread_4, kinetic_energy_data_SN_run_1_spread_4 = calendering_plot_processing(
        simulation_directory_SN_run_1_spread_4)


    # ==LOADING AND UNLOADING HEIGHTS===================================================================================

    loading_height_SN_run_1,unload_height_SN_run_1, specific_elastic_springback_SN_run_1 = load_unload_height_func(calendering_surface_pressure_SN_run_1,
                                                                             calendering_surface_position_SN_run_1)
    loading_height_SN_run_1_spread_1,unload_height_SN_run_1_spread_1, specific_elastic_springback_SN_run_1_spread_1 = load_unload_height_func(calendering_surface_pressure_SN_run_1_spread_1,
                                                                             calendering_surface_position_SN_run_1_spread_1)
    loading_height_SN_run_1_spread_2,unload_height_SN_run_1_spread_2, specific_elastic_springback_SN_run_1_spread_2 = load_unload_height_func(calendering_surface_pressure_SN_run_1_spread_2,
                                                                             calendering_surface_position_SN_run_1_spread_2)
    loading_height_SN_run_1_spread_3,unload_height_SN_run_1_spread_3, specific_elastic_springback_SN_run_1_spread_3 = load_unload_height_func(calendering_surface_pressure_SN_run_1_spread_3,
                                                                             calendering_surface_position_SN_run_1_spread_3)
    loading_height_SN_run_1_spread_4,unload_height_SN_run_1_spread_4, specific_elastic_springback_SN_run_1_spread_4 = load_unload_height_func(calendering_surface_pressure_SN_run_1_spread_4,
                                                                             calendering_surface_position_SN_run_1_spread_4)
    # print(loading_height_SN_run_1, unload_height_SN_run_1,elastic_recovery_SN_run_1)
    # print(loading_height_SN_run_1_spread_1,unload_height_SN_run_1_spread_1, elastic_recovery_SN_run_1_spread_1)
    # print(loading_height_SN_run_1_spread_2, unload_height_SN_run_1_spread_2, elastic_recovery_SN_run_1_spread_2)
    # print(loading_height_SN_run_1_spread_3, unload_height_SN_run_1_spread_3, elastic_recovery_SN_run_1_spread_3)
    # print(loading_height_SN_run_1_spread_4, unload_height_SN_run_1_spread_4, elastic_recovery_SN_run_1_spread_4)
    print('Mean of elastic recovery: ')
    print(np.mean([specific_elastic_springback_SN_run_1,specific_elastic_springback_SN_run_1_spread_1,specific_elastic_springback_SN_run_1_spread_2,specific_elastic_springback_SN_run_1_spread_3,specific_elastic_springback_SN_run_1_spread_4]))
    print('Standard deviation of elastic recovery: ')
    print(np.std([specific_elastic_springback_SN_run_1,specific_elastic_springback_SN_run_1_spread_1,specific_elastic_springback_SN_run_1_spread_2,specific_elastic_springback_SN_run_1_spread_3,specific_elastic_springback_SN_run_1_spread_4]))

    #==LOADING PLOT===============================================================================================

    max_initial_height = max([calendering_surface_position_SN_run_1[0],calendering_surface_position_SN_run_1_spread_1[0],calendering_surface_position_SN_run_1_spread_2[0],calendering_surface_position_SN_run_1_spread_3[0],calendering_surface_position_SN_run_1_spread_4[0]])
    max_unload_height = max([calendering_surface_position_SN_run_1[-1],calendering_surface_position_SN_run_1_spread_1[-1],calendering_surface_position_SN_run_1_spread_2[-1],calendering_surface_position_SN_run_1_spread_3[-1],calendering_surface_position_SN_run_1_spread_4[-1]])
    min_calendering_height = 1.05
    common_array_loading = np.linspace(min_calendering_height,max_initial_height,1000)
    common_array_unloading = np.linspace(min_calendering_height,max_unload_height)

    max_displacement_index_SN_run_1 = max_displacement_index_func(calendering_surface_position_SN_run_1)
    calendering_surface_position_loading_SN_run_1 = calendering_surface_position_SN_run_1[:max_displacement_index_SN_run_1]
    calendering_surface_pressure_loading_SN_run_1 = calendering_surface_pressure_SN_run_1[:max_displacement_index_SN_run_1]
    calendering_surface_pressure_loading_SN_run_1_int_pol = interpolate.interp1d(calendering_surface_position_loading_SN_run_1,calendering_surface_pressure_loading_SN_run_1,bounds_error=False,fill_value=0)
    calendering_surface_pressure_loading_SN_run_1_common = calendering_surface_pressure_loading_SN_run_1_int_pol(common_array_loading)

    calendering_surface_position_unloading_SN_run_1 = calendering_surface_position_SN_run_1[max_displacement_index_SN_run_1:]
    calendering_surface_pressure_unloading_SN_run_1 = calendering_surface_pressure_SN_run_1[max_displacement_index_SN_run_1:]
    calendering_surface_pressure_unloading_SN_run_1_int_pol = interpolate.interp1d(calendering_surface_position_unloading_SN_run_1,calendering_surface_pressure_unloading_SN_run_1,bounds_error=False,fill_value=0)
    calendering_surface_pressure_unloading_SN_run_1_common = calendering_surface_pressure_unloading_SN_run_1_int_pol(common_array_unloading)

    max_displacement_index_SN_run_1_spread_1 = max_displacement_index_func(calendering_surface_position_SN_run_1_spread_1)
    calendering_surface_position_loading_SN_run_1_spread_1 = calendering_surface_position_SN_run_1_spread_1[:max_displacement_index_SN_run_1_spread_1]
    calendering_surface_pressure_loading_SN_run_1_spread_1 = calendering_surface_pressure_SN_run_1_spread_1[:max_displacement_index_SN_run_1_spread_1]
    calendering_surface_pressure_loading_SN_run_1_spread_1_int_pol = interpolate.interp1d(
        calendering_surface_position_loading_SN_run_1_spread_1, calendering_surface_pressure_loading_SN_run_1_spread_1,bounds_error=False, fill_value=0)
    calendering_surface_pressure_loading_SN_run_1_spread_1_common = calendering_surface_pressure_loading_SN_run_1_spread_1_int_pol(
        common_array_loading)

    calendering_surface_position_unloading_SN_run_1_spread_1 = calendering_surface_position_SN_run_1_spread_1[
                                                      max_displacement_index_SN_run_1_spread_1:]
    calendering_surface_pressure_unloading_SN_run_1_spread_1 = calendering_surface_pressure_SN_run_1_spread_1[
                                                      max_displacement_index_SN_run_1_spread_1:]
    calendering_surface_pressure_unloading_SN_run_1_spread_1_int_pol = interpolate.interp1d(
        calendering_surface_position_unloading_SN_run_1_spread_1, calendering_surface_pressure_unloading_SN_run_1_spread_1,
        bounds_error=False, fill_value=0)
    calendering_surface_pressure_unloading_SN_run_1_spread_1_common = calendering_surface_pressure_unloading_SN_run_1_spread_1_int_pol(
        common_array_unloading)

    max_displacement_index_SN_run_1_spread_2 = max_displacement_index_func(calendering_surface_position_SN_run_1_spread_2)
    calendering_surface_position_loading_SN_run_1_spread_2 = calendering_surface_position_SN_run_1_spread_2[:max_displacement_index_SN_run_1_spread_2]
    calendering_surface_pressure_loading_SN_run_1_spread_2 = calendering_surface_pressure_SN_run_1_spread_2[:max_displacement_index_SN_run_1_spread_2]
    calendering_surface_pressure_loading_SN_run_1_spread_2_int_pol = interpolate.interp1d(
        calendering_surface_position_loading_SN_run_1_spread_2, calendering_surface_pressure_loading_SN_run_1_spread_2,bounds_error=False,fill_value=0)
    calendering_surface_pressure_loading_SN_run_1_spread_2_common = calendering_surface_pressure_loading_SN_run_1_spread_2_int_pol(
        common_array_loading)

    calendering_surface_position_unloading_SN_run_1_spread_2 = calendering_surface_position_SN_run_1_spread_2[
                                                      max_displacement_index_SN_run_1_spread_2:]
    calendering_surface_pressure_unloading_SN_run_1_spread_2 = calendering_surface_pressure_SN_run_1_spread_2[
                                                      max_displacement_index_SN_run_1_spread_2:]
    calendering_surface_pressure_unloading_SN_run_1_spread_2_int_pol = interpolate.interp1d(
        calendering_surface_position_unloading_SN_run_1_spread_2, calendering_surface_pressure_unloading_SN_run_1_spread_2,
        bounds_error=False, fill_value=0)
    calendering_surface_pressure_unloading_SN_run_1_spread_2_common = calendering_surface_pressure_unloading_SN_run_1_spread_2_int_pol(
        common_array_unloading)

    max_displacement_index_SN_run_1_spread_3 = max_displacement_index_func(calendering_surface_position_SN_run_1_spread_3)
    calendering_surface_position_loading_SN_run_1_spread_3 = calendering_surface_position_SN_run_1_spread_3[:max_displacement_index_SN_run_1_spread_3]
    calendering_surface_pressure_loading_SN_run_1_spread_3 = calendering_surface_pressure_SN_run_1_spread_3[:max_displacement_index_SN_run_1_spread_3]
    calendering_surface_pressure_loading_SN_run_1_spread_3_int_pol = interpolate.interp1d(
        calendering_surface_position_loading_SN_run_1_spread_3, calendering_surface_pressure_loading_SN_run_1_spread_3,bounds_error=False,fill_value=0)
    calendering_surface_pressure_loading_SN_run_1_spread_3_common = calendering_surface_pressure_loading_SN_run_1_spread_3_int_pol(
        common_array_loading)

    calendering_surface_position_unloading_SN_run_1_spread_3 = calendering_surface_position_SN_run_1_spread_3[
                                                      max_displacement_index_SN_run_1_spread_3:]
    calendering_surface_pressure_unloading_SN_run_1_spread_3 = calendering_surface_pressure_SN_run_1_spread_3[
                                                      max_displacement_index_SN_run_1_spread_3:]
    calendering_surface_pressure_unloading_SN_run_1_spread_3_int_pol = interpolate.interp1d(
        calendering_surface_position_unloading_SN_run_1_spread_3, calendering_surface_pressure_unloading_SN_run_1_spread_3,
        bounds_error=False, fill_value=0)
    calendering_surface_pressure_unloading_SN_run_1_spread_3_common = calendering_surface_pressure_unloading_SN_run_1_spread_3_int_pol(
        common_array_unloading)


    max_displacement_index_SN_run_1_spread_4 = max_displacement_index_func(calendering_surface_position_SN_run_1_spread_4)
    calendering_surface_position_loading_SN_run_1_spread_4 = calendering_surface_position_SN_run_1_spread_4[:max_displacement_index_SN_run_1_spread_4]
    calendering_surface_pressure_loading_SN_run_1_spread_4 = calendering_surface_pressure_SN_run_1_spread_4[:max_displacement_index_SN_run_1_spread_4]
    calendering_surface_pressure_loading_SN_run_1_spread_4_int_pol = interpolate.interp1d(
        calendering_surface_position_loading_SN_run_1_spread_4, calendering_surface_pressure_loading_SN_run_1_spread_4,bounds_error=False,fill_value=0)
    calendering_surface_pressure_loading_SN_run_1_spread_4_common = calendering_surface_pressure_loading_SN_run_1_spread_4_int_pol(
        common_array_loading)

    calendering_surface_position_unloading_SN_run_1_spread_4 = calendering_surface_position_SN_run_1_spread_4[
                                                      max_displacement_index_SN_run_1_spread_4:]
    calendering_surface_pressure_unloading_SN_run_1_spread_4 = calendering_surface_pressure_SN_run_1_spread_4[
                                                      max_displacement_index_SN_run_1_spread_4:]
    calendering_surface_pressure_unloading_SN_run_1_spread_4_int_pol = interpolate.interp1d(
        calendering_surface_position_unloading_SN_run_1_spread_4, calendering_surface_pressure_unloading_SN_run_1_spread_4,
        bounds_error=False, fill_value=0)
    calendering_surface_pressure_unloading_SN_run_1_spread_4_common = calendering_surface_pressure_unloading_SN_run_1_spread_4_int_pol(
        common_array_unloading)

    all_pressures_loading = np.vstack((calendering_surface_pressure_loading_SN_run_1_common,
                                       calendering_surface_pressure_loading_SN_run_1_spread_1_common,
                                       calendering_surface_pressure_loading_SN_run_1_spread_2_common,
                                       calendering_surface_pressure_loading_SN_run_1_spread_3_common,
                                       calendering_surface_pressure_loading_SN_run_1_spread_4_common))
    all_pressures_unloading = np.vstack((calendering_surface_pressure_unloading_SN_run_1_common,
                                       calendering_surface_pressure_unloading_SN_run_1_spread_1_common,
                                       calendering_surface_pressure_unloading_SN_run_1_spread_2_common,
                                       calendering_surface_pressure_unloading_SN_run_1_spread_3_common,
                                       calendering_surface_pressure_unloading_SN_run_1_spread_4_common))

    std_loading = np.std(all_pressures_loading,axis=0)
    mean_loading = np.mean(all_pressures_loading,axis=0)

    std_unloading = np.std(all_pressures_unloading,axis=0)
    mean_unloading = np.mean(all_pressures_unloading,axis=0)

    fig_calendering_loading, ax_calendering_loading = plt.subplots()
#    ax_calendering_loading.plot(calendering_surface_position_SN_run_1[:max_displacement_index_SN_run_1 ] * 1E2,
#                                 calendering_surface_pressure_SN_run_1[:max_displacement_index_SN_run_1 ] * 1E-6, label=r'$Run 1$')
#
#    ax_calendering_loading.plot(calendering_surface_position_SN_run_1[max_displacement_index_SN_run_1:] * 1E2,
#                                calendering_surface_pressure_SN_run_1[max_displacement_index_SN_run_1:] * 1E-6,
#    ax_calendering_loading.plot(calendering_surface_position_SN_run_1_spread_1[max_displacement_index_SN_run_1_spread_1:] * 1E2,
#                                calendering_surface_pressure_SN_run_1_spread_1[max_displacement_index_SN_run_1_spread_1:] * 1E-6,
#                                label=r'$Run 2$')
#    ax_calendering_loading.plot(calendering_surface_position_SN_run_1_spread_2[max_displacement_index_SN_run_1_spread_2:] * 1E2,
#                                calendering_surface_pressure_SN_run_1_spread_2[max_displacement_index_SN_run_1_spread_2:] * 1E-6,
#                                label=r'$Run 3$')


#    ax_calendering_loading.plot(calendering_surface_position_SN_run_1_spread_1[:max_displacement_index_SN_run_1_spread_1] * 1E2,
#                                 calendering_surface_pressure_SN_run_1_spread_1[:max_displacement_index_SN_run_1_spread_1] * 1E-6, label=r'$Run 2$')

#    ax_calendering_loading.plot(calendering_surface_position_SN_run_1_spread_2[:max_displacement_index_SN_run_1_spread_2] * 1E2,
#                                 calendering_surface_pressure_SN_run_1_spread_2[:max_displacement_index_SN_run_1_spread_2] * 1E-6, label=r'$Run 3$')

    lns_mean_loading = ax_calendering_loading.plot(common_array_loading * 1E2, mean_loading * 1E-6,'k-')
    lns_std_loading = plt.fill_between(common_array_loading * 1E2, (mean_loading-std_loading) * 1E-6,
                                       (mean_loading + std_loading) * 1E-6, lw=2, color='C0', alpha=1)

    ax_calendering_loading.plot(common_array_unloading * 1E2, mean_unloading * 1E-6, 'k-')
    plt.fill_between(common_array_unloading * 1E2, (mean_unloading - std_unloading) * 1E-6,
                     (mean_unloading + std_unloading) * 1E-6, color='C0',alpha=1)

    ax_calendering_loading.set_ylim(ymin=0)
    ax_calendering_loading.set_xlim(xmin=104.8,xmax=125)
    ax_calendering_loading.xaxis.set_major_locator(MultipleLocator(5))

    ax_calendering_loading.legend([lns_mean_loading[0],lns_std_loading],['Mean','Standard deviation'],loc='best')
    ax_calendering_loading.set_ylabel("Calendering surface pressure [MPa]")
    ax_calendering_loading.set_xlabel("Calendering surface height [µm]")

    fname = fig_dir + 'calendering_surface_pressure_surface_position_full_standard_deviation'
    plt.savefig(fname)


    #==STATISTICAL SPREAD===============================================================================================
    fig_full_calendering_sim, ax_full_calendering_sim = plt.subplots()
    ax_full_calendering_sim.plot(calendering_surface_position_SN_run_1[:] * 1E2, calendering_surface_pressure_SN_run_1[:] * 1E-6, label=r'$Run 1$')
    ax_full_calendering_sim.plot(calendering_surface_position_SN_run_1_spread_1[:] * 1E2,
                                 calendering_surface_pressure_SN_run_1_spread_1[:] * 1E-6, label=r'$Run 2$')
    ax_full_calendering_sim.plot(calendering_surface_position_SN_run_1_spread_2[:] * 1E2,
                                 calendering_surface_pressure_SN_run_1_spread_2[:] * 1E-6, label=r'$Run 3$')
    ax_full_calendering_sim.plot(calendering_surface_position_SN_run_1_spread_3[:] * 1E2,
                                 calendering_surface_pressure_SN_run_1_spread_3[:] * 1E-6, label=r'$Run 4$')
    ax_full_calendering_sim.plot(calendering_surface_position_SN_run_1_spread_4[:] * 1E2,
                                 calendering_surface_pressure_SN_run_1_spread_4[:] * 1E-6, label=r'$Run 5$')
    ax_full_calendering_sim.set_xlim(xmin=104.8,xmax=125)
    ax_full_calendering_sim.xaxis.set_major_locator(MultipleLocator(5))

    ax_full_calendering_sim.set_ylim(ymin=0)
    ax_full_calendering_sim.set_ylabel("Calendering surface pressure [MPa]")
    ax_full_calendering_sim.set_xlabel("Calendering surface height [µm]")
    #    ax_full_calendering_sim.set_title('calendering surface pressure')
    fig_full_calendering_sim.tight_layout()
    ax_full_calendering_sim.legend(loc="best")

    fname = fig_dir + 'calendering_surface_pressure_surface_position_full_statistical_spread'
    plt.savefig(fname)

    plt.show()
