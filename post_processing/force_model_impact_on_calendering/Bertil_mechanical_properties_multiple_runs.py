from force_model_impact_on_calendering.Bertil_calendering_pressure import bertil_data_gatherer, one_file_reader, contact_counter_bertil

import numpy as np
import matplotlib.pyplot as plt
import shutil
import os
import matplotlib
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import pandas as pd


def unload_index_locator(linear_strains):
    break_index = np.where(linear_strains == np.abs(linear_strains).max())[-1]

    return break_index


def time_duplicate_remover(time_list, linear_strain_list):
    rm_list = []
    for i in range(len(time_list) - 1):
        if time_list[i] == time_list[i + 1]:
            rm_list.append(i)
    time_short = np.delete(time_list, rm_list)
    linear_strain_short = np.delete(linear_strain_list, rm_list)
    return time_short, linear_strain_short


def stress_and_linear_strain_finder(periodic_BC_data, force_fabric_tensor_data, surface_position_data):
    time = periodic_BC_data[:, 0]

    x_side_length = periodic_BC_data[:, 2] - periodic_BC_data[:, 1]
    x_side_length_0 = x_side_length[0]
    y_side_length = periodic_BC_data[:, 4] - periodic_BC_data[:, 3]

    # t0 = 1.11  # Chould be chaged later, defined in another way?======================================================
    t0 = float(surface_position_data[0, 32]) - 1
    vol = x_side_length * y_side_length * t0
    sxx = -force_fabric_tensor_data[:, 1] / vol
    syy = -force_fabric_tensor_data[:, 5] / vol
    szz = -force_fabric_tensor_data[:, 9] / vol

    tau_xy = -force_fabric_tensor_data[:, 2] / vol
    tau_xz = -force_fabric_tensor_data[:, 3] / vol

    tau_yx = -force_fabric_tensor_data[:, 4] / vol
    tau_yz = -force_fabric_tensor_data[:, 6] / vol

    tau_zx = -force_fabric_tensor_data[:, 7] / vol
    tau_zy = -force_fabric_tensor_data[:, 8] / vol

    linear_strain = (x_side_length[:] - x_side_length_0) / x_side_length_0

    return time, linear_strain, sxx, syy, szz, tau_xy, tau_xz, tau_yx, tau_yz, tau_zx, tau_zy


def stiffness_finder(stress_vec, strain_vec):
    count = 0
    n = len(strain_vec)
    max_index = []
    min_index = []
    # Tension
    if strain_vec[1] > 0:
        for i in range(1, n - 1):
            # maximum values
            if (strain_vec[i] > strain_vec[i + 1] and strain_vec[i] >= strain_vec[i - 1]):
                max_index.append(i)
                count += 1
            # Minimum values
            if (strain_vec[i] <= strain_vec[i + 1] and strain_vec[i] < strain_vec[i - 1]):
                min_index.append(i)
                count += 1
            if count >= 10:
                break
    # Compression
    if strain_vec[1] < 0:
        for i in range(1, n - 1):
            # Maximum values
            if (strain_vec[i] >= strain_vec[i + 1] and strain_vec[i] > strain_vec[i - 1]):
                max_index.append(i)
                count += 1
            # Minimum values
            if (strain_vec[i] < strain_vec[i + 1] and strain_vec[i] <= strain_vec[i - 1]):
                min_index.append(i)
                count += 1
            if count >= 10:
                break
    stiffness_values = np.array([])
    for i in range(5):
        # ==================STRAIN VÄRDEN LÄGGS INTE TILL I ARRAYN...===========
        stiffness_values = np.append(stiffness_values, (stress_vec[max_index[i]] - stress_vec[min_index[i]]) / (
                    strain_vec[max_index[i]] - strain_vec[min_index[i]]))
    if strain_vec[1] < 0:
        strain_points = strain_vec[min_index]
    if strain_vec[1] > 0:
        strain_points = strain_vec[max_index]

    return strain_points, stiffness_values  # returns strain in [-] and stiffness in Pa


def mech_plot_prop(simulation_directory_SN_run_1):
    force_data_compression, surface_force_index_compression, surface_position_index_compression, \
    surface_position_data_compression, periodic_BC_data_compression, force_fabric_tensor_data_compression, \
    kinetic_energy_data_compression = bertil_data_gatherer(
        simulation_directory_SN_run_1 + '_compression')

    force_data_tension, surface_force_index_tension, surface_position_index_tension, surface_position_data_tension, \
    periodic_BC_data_tension, force_fabric_tensor_data_tension, kinetic_energy_data_tension = bertil_data_gatherer(
        simulation_directory_SN_run_1 + '_tension')

    time_tension, linear_strain_tension, sxx_tension, syy_tension, szz_tension, tau_xy_tension, tau_xz_tension, \
    tau_yx_tension, tau_yz_tension, tau_zx_tension, tau_zy_tension = \
        stress_and_linear_strain_finder(periodic_BC_data_tension,
                                        force_fabric_tensor_data_tension, surface_position_data_tension)
    time_compression, linear_strain_compression, sxx_compression, syy_compression, szz_compression, tau_xy_compression, \
    tau_xz_compression, tau_yx_compression, tau_yz_compression, tau_zx_compression, tau_zy_compression = \
        stress_and_linear_strain_finder(periodic_BC_data_compression, force_fabric_tensor_data_compression,
                                        surface_position_data_compression)

    return time_tension, linear_strain_tension, sxx_tension, syy_tension, szz_tension, tau_xy_tension, tau_xz_tension, \
           tau_yx_tension, tau_yz_tension, tau_zx_tension, tau_zy_tension, time_compression, linear_strain_compression, \
           sxx_compression, syy_compression, szz_compression, tau_xy_compression, \
           tau_xz_compression, tau_yx_compression, tau_yz_compression, tau_zx_compression, tau_zy_compression


def stiffness_func(sxx_tension, linear_strain_tension, sxx_compression, linear_strain_compression):
    strain_points_tension, stiffness_values_tension = stiffness_finder(sxx_tension, linear_strain_tension)
    strain_points_compression, stiffness_values_compression = stiffness_finder(sxx_compression,
                                                                               linear_strain_compression)

    strain_points_compression = np.flip(strain_points_compression)
    stiffness_values_compression = np.flip(stiffness_values_compression)

    strain_points_total = np.concatenate((strain_points_compression, strain_points_tension))
    stiffness_values_total = np.concatenate((stiffness_values_compression, stiffness_values_tension))
    return strain_points_total, stiffness_values_total


if __name__ == '__main__':

    simulation_directory_SN_run_1 = '/scratch/users/axlun/DEMsim/results/final_runs/SN_run_1/electrode_mechanical_loading_hertz'
    simulation_directory_SN_run_1_2E = '/scratch/users/axlun/DEMsim/results/final_runs/SN_run_1_2E/electrode_mechanical_loading_hertz'
    simulation_directory_SN_run_1_El_Pl = '/scratch/users/axlun/DEMsim/results/final_runs/SN_run_1_El_Pl/electrode_mechanical_loading_hertz'
    simulation_directory_SN_run_1_br_05 = '/scratch/users/axlun/DEMsim/results/final_runs/SN_run_1_br_05/electrode_mechanical_loading_hertz'
    simulation_directory_SN_run_1_br_10 = '/scratch/users/axlun/DEMsim/results/final_runs/SN_run_1_br_10/electrode_mechanical_loading_hertz'
    simulation_directory_SN_run_1_br_175 = '/scratch/users/axlun/DEMsim/results/final_runs/SN_run_1_br_175/electrode_mechanical_loading_hertz'
    simulation_directory_SN_run_1_EL_particle = '/scratch/users/axlun/DEMsim/results/final_runs/SN_run_1_EL_particle/electrode_mechanical_loading_hertz'
    simulation_directory_SN_run_1_q_0_rmin_3_rmax_10 = '/scratch/users/axlun/DEMsim/results/final_runs/SN_run_1_q_0_rmin_3_rmax_10/electrode_mechanical_loading_hertz'
    simulation_directory_SN_run_1_q_0_rmin_4_rmax_10 = '/scratch/users/axlun/DEMsim/results/final_runs/SN_run_1_q_0_rmin_4_rmax_10/electrode_mechanical_loading_hertz'

    time_tension_SN_run_1, linear_strain_tension_SN_run_1, sxx_tension_SN_run_1, syy_tension_SN_run_1, szz_tension_SN_run_1, tau_xy_tension_SN_run_1, tau_xz_tension_SN_run_1, \
    tau_yx_tension_SN_run_1, tau_yz_tension_SN_run_1, tau_zx_tension_SN_run_1, tau_zy_tension_SN_run_1, time_compression_SN_run_1, linear_strain_compression_SN_run_1, \
    sxx_compression_SN_run_1, syy_compression_SN_run_1, szz_compression_SN_run_1, tau_xy_compression_SN_run_1, \
    tau_xz_compression_SN_run_1, tau_yx_compression_SN_run_1, tau_yz_compression_SN_run_1, tau_zx_compression_SN_run_1, tau_zy_compression_SN_run_1 = mech_plot_prop(simulation_directory_SN_run_1)

    time_tension_SN_run_1_2E, linear_strain_tension_SN_run_1_2E, sxx_tension_SN_run_1_2E, syy_tension_SN_run_1_2E, szz_tension_SN_run_1_2E, tau_xy_tension_SN_run_1_2E, tau_xz_tension_SN_run_1_2E, \
    tau_yx_tension_SN_run_1_2E, tau_yz_tension_SN_run_1_2E, tau_zx_tension_SN_run_1_2E, tau_zy_tension_SN_run_1_2E, time_compression_SN_run_1_2E, linear_strain_compression_SN_run_1_2E, \
    sxx_compression_SN_run_1_2E, syy_compression_SN_run_1_2E, szz_compression_SN_run_1_2E, tau_xy_compression_SN_run_1_2E, \
    tau_xz_compression_SN_run_1_2E, tau_yx_compression_SN_run_1_2E, tau_yz_compression_SN_run_1_2E, tau_zx_compression_SN_run_1_2E, tau_zy_compression_SN_run_1_2E = mech_plot_prop(
        simulation_directory_SN_run_1_2E)

    time_tension_SN_run_1_El_Pl, linear_strain_tension_SN_run_1_El_Pl, sxx_tension_SN_run_1_El_Pl, syy_tension_SN_run_1_El_Pl, szz_tension_SN_run_1_El_Pl, tau_xy_tension_SN_run_1_El_Pl, tau_xz_tension_SN_run_1_El_Pl, \
    tau_yx_tension_SN_run_1_El_Pl, tau_yz_tension_SN_run_1_El_Pl, tau_zx_tension_SN_run_1_El_Pl, tau_zy_tension_SN_run_1_El_Pl, time_compression_SN_run_1_El_Pl, linear_strain_compression_SN_run_1_El_Pl, \
    sxx_compression_SN_run_1_El_Pl, syy_compression_SN_run_1_El_Pl, szz_compression_SN_run_1_El_Pl, tau_xy_compression_SN_run_1_El_Pl, \
    tau_xz_compression_SN_run_1_El_Pl, tau_yx_compression_SN_run_1_El_Pl, tau_yz_compression_SN_run_1_El_Pl, tau_zx_compression_SN_run_1_El_Pl, tau_zy_compression_SN_run_1_El_Pl = mech_plot_prop(
        simulation_directory_SN_run_1_El_Pl)

    time_tension_SN_run_1_br_05, linear_strain_tension_SN_run_1_br_05, sxx_tension_SN_run_1_br_05, syy_tension_SN_run_1_br_05, szz_tension_SN_run_1_br_05, tau_xy_tension_SN_run_1_br_05, tau_xz_tension_SN_run_1_br_05, \
    tau_yx_tension_SN_run_1_br_05, tau_yz_tension_SN_run_1_br_05, tau_zx_tension_SN_run_1_br_05, tau_zy_tension_SN_run_1_br_05, time_compression_SN_run_1_br_05, linear_strain_compression_SN_run_1_br_05, \
    sxx_compression_SN_run_1_br_05, syy_compression_SN_run_1_br_05, szz_compression_SN_run_1_br_05, tau_xy_compression_SN_run_1_br_05, \
    tau_xz_compression_SN_run_1_br_05, tau_yx_compression_SN_run_1_br_05, tau_yz_compression_SN_run_1_br_05, tau_zx_compression_SN_run_1_br_05, tau_zy_compression_SN_run_1_br_05 = mech_plot_prop(
        simulation_directory_SN_run_1_br_05)

    time_tension_SN_run_1_br_10, linear_strain_tension_SN_run_1_br_10, sxx_tension_SN_run_1_br_10, syy_tension_SN_run_1_br_10, szz_tension_SN_run_1_br_10, tau_xy_tension_SN_run_1_br_10, tau_xz_tension_SN_run_1_br_10, \
    tau_yx_tension_SN_run_1_br_10, tau_yz_tension_SN_run_1_br_10, tau_zx_tension_SN_run_1_br_10, tau_zy_tension_SN_run_1_br_10, time_compression_SN_run_1_br_10, linear_strain_compression_SN_run_1_br_10, \
    sxx_compression_SN_run_1_br_10, syy_compression_SN_run_1_br_10, szz_compression_SN_run_1_br_10, tau_xy_compression_SN_run_1_br_10, \
    tau_xz_compression_SN_run_1_br_10, tau_yx_compression_SN_run_1_br_10, tau_yz_compression_SN_run_1_br_10, tau_zx_compression_SN_run_1_br_10, tau_zy_compression_SN_run_1_br_10 = mech_plot_prop(
        simulation_directory_SN_run_1_br_10)

    time_tension_SN_run_1_br_175, linear_strain_tension_SN_run_1_br_175, sxx_tension_SN_run_1_br_175, syy_tension_SN_run_1_br_175, szz_tension_SN_run_1_br_175, tau_xy_tension_SN_run_1_br_175, tau_xz_tension_SN_run_1_br_175, \
    tau_yx_tension_SN_run_1_br_175, tau_yz_tension_SN_run_1_br_175, tau_zx_tension_SN_run_1_br_175, tau_zy_tension_SN_run_1_br_175, time_compression_SN_run_1_br_175, linear_strain_compression_SN_run_1_br_175, \
    sxx_compression_SN_run_1_br_175, syy_compression_SN_run_1_br_175, szz_compression_SN_run_1_br_175, tau_xy_compression_SN_run_1_br_175, \
    tau_xz_compression_SN_run_1_br_175, tau_yx_compression_SN_run_1_br_175, tau_yz_compression_SN_run_1_br_175, tau_zx_compression_SN_run_1_br_175, tau_zy_compression_SN_run_1_br_175 = mech_plot_prop(
        simulation_directory_SN_run_1_br_175)

    time_tension_SN_run_1_EL_particle, linear_strain_tension_SN_run_1_EL_particle, sxx_tension_SN_run_1_EL_particle, syy_tension_SN_run_1_EL_particle, szz_tension_SN_run_1_EL_particle, tau_xy_tension_SN_run_1_EL_particle, tau_xz_tension_SN_run_1_EL_particle, \
    tau_yx_tension_SN_run_1_EL_particle, tau_yz_tension_SN_run_1_EL_particle, tau_zx_tension_SN_run_1_EL_particle, tau_zy_tension_SN_run_1_EL_particle, time_compression_SN_run_1_EL_particle, linear_strain_compression_SN_run_1_EL_particle, \
    sxx_compression_SN_run_1_EL_particle, syy_compression_SN_run_1_EL_particle, szz_compression_SN_run_1_EL_particle, tau_xy_compression_SN_run_1_EL_particle, \
    tau_xz_compression_SN_run_1_EL_particle, tau_yx_compression_SN_run_1_EL_particle, tau_yz_compression_SN_run_1_EL_particle, tau_zx_compression_SN_run_1_EL_particle, tau_zy_compression_SN_run_1_EL_particle = mech_plot_prop(
        simulation_directory_SN_run_1_EL_particle)

    time_tension_SN_run_1_q_0_rmin_3_rmax_10, linear_strain_tension_SN_run_1_q_0_rmin_3_rmax_10, sxx_tension_SN_run_1_q_0_rmin_3_rmax_10, syy_tension_SN_run_1_q_0_rmin_3_rmax_10, szz_tension_SN_run_1_q_0_rmin_3_rmax_10, tau_xy_tension_SN_run_1_q_0_rmin_3_rmax_10, tau_xz_tension_SN_run_1_q_0_rmin_3_rmax_10, \
    tau_yx_tension_SN_run_1_q_0_rmin_3_rmax_10, tau_yz_tension_SN_run_1_q_0_rmin_3_rmax_10, tau_zx_tension_SN_run_1_q_0_rmin_3_rmax_10, tau_zy_tension_SN_run_1_q_0_rmin_3_rmax_10, time_compression_SN_run_1_q_0_rmin_3_rmax_10, linear_strain_compression_SN_run_1_q_0_rmin_3_rmax_10, \
    sxx_compression_SN_run_1_q_0_rmin_3_rmax_10, syy_compression_SN_run_1_q_0_rmin_3_rmax_10, szz_compression_SN_run_1_q_0_rmin_3_rmax_10, tau_xy_compression_SN_run_1_q_0_rmin_3_rmax_10, \
    tau_xz_compression_SN_run_1_q_0_rmin_3_rmax_10, tau_yx_compression_SN_run_1_q_0_rmin_3_rmax_10, tau_yz_compression_SN_run_1_q_0_rmin_3_rmax_10, tau_zx_compression_SN_run_1_q_0_rmin_3_rmax_10, tau_zy_compression_SN_run_1_q_0_rmin_3_rmax_10 = mech_plot_prop(
        simulation_directory_SN_run_1_q_0_rmin_3_rmax_10)

    time_tension_SN_run_1_q_0_rmin_4_rmax_10, linear_strain_tension_SN_run_1_q_0_rmin_4_rmax_10, sxx_tension_SN_run_1_q_0_rmin_4_rmax_10, syy_tension_SN_run_1_q_0_rmin_4_rmax_10, szz_tension_SN_run_1_q_0_rmin_4_rmax_10, tau_xy_tension_SN_run_1_q_0_rmin_4_rmax_10, tau_xz_tension_SN_run_1_q_0_rmin_4_rmax_10, \
    tau_yx_tension_SN_run_1_q_0_rmin_4_rmax_10, tau_yz_tension_SN_run_1_q_0_rmin_4_rmax_10, tau_zx_tension_SN_run_1_q_0_rmin_4_rmax_10, tau_zy_tension_SN_run_1_q_0_rmin_4_rmax_10, time_compression_SN_run_1_q_0_rmin_4_rmax_10, linear_strain_compression_SN_run_1_q_0_rmin_4_rmax_10, \
    sxx_compression_SN_run_1_q_0_rmin_4_rmax_10, syy_compression_SN_run_1_q_0_rmin_4_rmax_10, szz_compression_SN_run_1_q_0_rmin_4_rmax_10, tau_xy_compression_SN_run_1_q_0_rmin_4_rmax_10, \
    tau_xz_compression_SN_run_1_q_0_rmin_4_rmax_10, tau_yx_compression_SN_run_1_q_0_rmin_4_rmax_10, tau_yz_compression_SN_run_1_q_0_rmin_4_rmax_10, tau_zx_compression_SN_run_1_q_0_rmin_4_rmax_10, tau_zy_compression_SN_run_1_q_0_rmin_4_rmax_10 = mech_plot_prop(
        simulation_directory_SN_run_1_q_0_rmin_4_rmax_10)

    stiffness_at_points_flag = 1

    fig_dir = 'C:/temp/figures/Bertil_mechanical_properties_multiple_runs/'
    try:
        shutil.rmtree(fig_dir)
        os.mkdir(fig_dir)
    except FileNotFoundError:
        print('Directory does not exist')
        os.mkdir(fig_dir)
    except:
        print('Directory could not be removed')
        quit()
    # ==PLOT PARAMETERS=====================================================================================================
    plt.style.use('axel_style')

    # ==EXPERIMENTAL DATA===================================================================================================
    exp_strain_points_compression = [-1.00, -1.10, -1.23, -1.41, -1.65]
    """
    #Wrong values read from from experiment result table
    Modulus_eps_dot_01_compression = [1.20, 1.04, 0.96, 1.19, 1.46]
    Modulus_eps_dot_05_compression = [1.43, 0.86, 0.97, 1.39, 1.79]
    Modulus_eps_dot_10_compression = [1.51, 1.79, 1.99, 2.09, 2.74]
    Modulus_eps_dot_100_compression = [1.55, 1.72, 1.70, 1.81, 2.28]
    Modulus_eps_dot_300_compression = [1.99, 1.52, 2.21, 2.80, 2.33]
    """
    Modulus_eps_dot_01_compression = [1.20, 1.43, 1.51, 1.55, 1.99]
    Modulus_eps_dot_05_compression = [1.04, 0.86, 1.79, 1.72, 1.52]
    Modulus_eps_dot_10_compression = [0.96, 0.97, 1.99, 1.70, 2.21]
    Modulus_eps_dot_100_compression = [1.19, 1.39, 2.09, 1.81, 2.80]
    Modulus_eps_dot_300_compression = [1.46, 1.79, 2.74, 2.28, 2.33]


    exp_strain_points_tension = [1.00, 1.10, 1.23, 1.41, 1.65]
    """
    #Wrong values read from from experiment result table
    Modulus_eps_dot_01_tension = [0.95, 0.50, 0.89, 0.61, 0.57]
    Modulus_eps_dot_05_tension = [0.78, 0.39, 1.17, 0.55, 0.46]
    Modulus_eps_dot_10_tension = [0.84, 0.59, 1.09, 0.59, 0.52]
    Modulus_eps_dot_100_tension = [0.76, 0.55, 0.91, 0.63, 0.67]
    Modulus_eps_dot_300_tension = [1.06, 0.55, 0.90, 0.95, 0.85]
    """
    Modulus_eps_dot_01_tension = [0.95, 0.78, 0.84, 0.76, 1.06]
    Modulus_eps_dot_05_tension = [0.50, 0.39, 0.59, 0.55, 0.55]
    Modulus_eps_dot_10_tension = [0.89, 1.17, 1.09, 0.91, 0.90]
    Modulus_eps_dot_100_tension = [0.61, 0.55, 0.59, 0.63, 0.95]
    Modulus_eps_dot_300_tension = [0.57, 0.46, 0.52, 0.67, 0.85]

    exp_strain_points = exp_strain_points_compression[::-1] + exp_strain_points_tension
    Modulus_eps_dot_01 = Modulus_eps_dot_01_compression[::-1] + Modulus_eps_dot_01_tension
    Modulus_eps_dot_05 = Modulus_eps_dot_05_compression[::-1] + Modulus_eps_dot_05_tension
    Modulus_eps_dot_10 = Modulus_eps_dot_10_compression[::-1] + Modulus_eps_dot_10_tension
    Modulus_eps_dot_100 = Modulus_eps_dot_100_compression[::-1] + Modulus_eps_dot_100_tension
    Modulus_eps_dot_300 = Modulus_eps_dot_300_compression[::-1] + Modulus_eps_dot_300_tension

    # # ==FIG 1 STRESS STRAIN IN COMPRESSION==============================================================================
    # fig_compression,ax_compression = plt.subplots()
    # ax_compression.set_ylabel('Stress [MPa]')
    # ax_compression.set_xlabel('Strain [%]')
    # lns_compression_xx = ax_compression.plot(linear_strain_compression[:] * 100, sxx_compression[:] / 1e6,
    #                                          label=r'$\sigma_{xx}$')
    # lns_compression_yy = ax_compression.plot(linear_strain_compression[:] * 100, syy_compression[:] / 1e6,
    #                                          label=r'$\sigma_{yy}$')
    # lns_compression_zz = ax_compression.plot(linear_strain_compression[:] * 100, szz_compression[:] / 1e6,
    #                                          label=r'$\sigma_{zz}$')
    # lns_compression_tau_xy = ax_compression.plot(linear_strain_compression[:] * 100, tau_xy_compression[:] / 1e6, label=r'$\tau_{xy}$')
    # lns_compression_tau_xz = ax_compression.plot(linear_strain_compression[:] * 100, tau_xz_compression[:] / 1e6, label=r'$\tau_{xz}$')
    # lns_compression_tau_yx = ax_compression.plot(linear_strain_compression[:] * 100, tau_yx_compression[:] / 1e6, label=r'$\tau_{yx}$')
    # lns_compression_tau_yz = ax_compression.plot(linear_strain_compression[:] * 100, tau_yz_compression[:] / 1e6, label=r'$\tau_{yz}$')
    # lns_compression_tau_zx = ax_compression.plot(linear_strain_compression[:] * 100, tau_zx_compression[:] / 1e6, label=r'$\tau_{zx}$')
    # lns_compression_tau_zy = ax_compression.plot(linear_strain_compression[:] * 100, tau_zy_compression[:] / 1e6, label=r'$\tau_{zy}$')
    #
    #
    # ax_compression.set_title('Macroscopic stress in compression')
    # plt.legend(loc='best')
    #
    # fname = fig_dir + 'compression_stress_strain'
    # plt.savefig(fname)
    #
    # # ==FIG 2 STRESS STRAIN IN TENSION==================================================================================
    # fig_tension,ax_tension = plt.subplots()
    # ax_tension.set_ylabel('Stress [MPa]')
    # ax_tension.set_xlabel('Strain [%]')
    # lns_tension_xx = ax_tension.plot(linear_strain_tension[:] * 100, sxx_tension[:] / 1e6, label=r'$\sigma_{xx}$')
    # lns__tension_yy = ax_tension.plot(linear_strain_tension[:] * 100, syy_tension[:] / 1e6, label=r'$\sigma_{yy}$')
    # lns_tension_zz = ax_tension.plot(linear_strain_tension[:] * 100, szz_tension[:] / 1e6, label=r'$\sigma_{zz}$')
    #
    # lns_tension_tau_xy = ax_tension.plot(linear_strain_tension[:] * 100, tau_xy_tension[:] / 1e6, label=r'$\tau_{xy}$')
    # lns_tension_tau_xz = ax_tension.plot(linear_strain_tension[:] * 100, tau_xz_tension[:] / 1e6, label=r'$\tau_{xz}$')
    # lns_tension_tau_yx = ax_tension.plot(linear_strain_tension[:] * 100, tau_yx_tension[:] / 1e6, label=r'$\tau_{yx}$')
    # lns_tension_tau_yz = ax_tension.plot(linear_strain_tension[:] * 100, tau_yz_tension[:] / 1e6, label=r'$\tau_{yz}$')
    # lns_tension_tau_zx = ax_tension.plot(linear_strain_tension[:] * 100, tau_zx_tension[:] / 1e6, label=r'$\tau_{zx}$')
    # lns_tension_tau_zy = ax_tension.plot(linear_strain_tension[:] * 100, tau_zy_tension[:] / 1e6, label=r'$\tau_{zy}$')
    #
    #
    # ax_tension.set_title('Macroscopic stress in tension')
    # plt.legend(loc='best')
    #
    # fname = fig_dir + 'tension_stress_strain'
    # plt.savefig(fname)
    # # ==FIG 3 STRESS AND STRAIN TO TIME FOR TENSION=====================================================================
    # fig_tension_time, ax_tension_time = plt.subplots()
    # ax_tension_time.set_ylabel("Stress [MPa]")
    # ax_tension_time.set_xlabel("time [s]")
    # lns_tension_xx_time = ax_tension_time.plot(time_tension, sxx_tension[:] / 1e6, label=r'$\sigma_{xx}$')
    # lns_tension_yy_time = ax_tension_time.plot(time_tension, syy_tension[:] / 1e6,
    #                                                    label=r'$\sigma_{yy}$')
    # lns_tension_zz_time = ax_tension_time.plot(time_tension, szz_tension[:] / 1e6,
    #                                                    label=r'$\sigma_{zz}$')
    # ax_tension_time_2 = ax_tension_time.twinx()
    # ##    ax2.set_xlabel("Calendering surface position [m]")
    # #    ax2.set_ylabel("Calendering surface position [m]")
    # lns_tension_time_2 = ax_tension_time_2.plot(time_tension, linear_strain_tension * 100, 'c', label=r'$\varepsilon_{xx}$')
    # ax_tension_time_2.set_ylabel("Strain [%]")
    # # lns3 = ax2.plot(time[max_index], linear_strain[max_index],'b',label='Position',linestyle="None",marker='x',markeredgewidth = 3)
    # # lns4 = ax2.plot(time[min_index], linear_strain[min_index], 'b', label='Position', linestyle="None", marker='x',markeredgewidth = 3)
    #
    # lns_tension = lns_tension_xx_time + lns_tension_yy_time+ lns_tension_zz_time + lns_tension_time_2
    # labs_tension = [l.get_label() for l in lns_tension]
    # ax_tension_time.legend(lns_tension, labs_tension, loc=0)
    # ax_tension_time.set_title('Macroscopic stress and strain in tension')
    #
    # fname = fig_dir + 'tension_stress_strain_time'
    # plt.savefig(fname)
    #
    # # ==FIG 4 STRESS AND STRAIN TO TIME FOR COMPRESSION=================================================================
    # fig_compression_time, ax_compression_time = plt.subplots()
    # ax_compression_time.set_ylabel("Stress [MPa]")
    # ax_compression_time.set_xlabel("time [s]")
    # lns_compression_xx_time = ax_compression_time.plot(time_compression, sxx_compression[:] / 1e6, label=r'$\sigma_{xx}$')
    # lns_compression_yy_time = ax_compression_time.plot(time_compression, syy_compression[:] / 1e6,
    #                                                    label=r'$\sigma_{yy}$')
    # lns_compression_zz_time = ax_compression_time.plot(time_compression, szz_compression[:] / 1e6, label=r'$\sigma_{zz}$')
    # ax_compression_time_2 = ax_compression_time.twinx()
    # lns_compression_time_2 = ax_compression_time_2.plot(time_compression, linear_strain_compression * 100, 'c', label=r'$\varepsilon_{xx}$')
    # ax_compression_time_2.set_ylabel("Strain [%]")
    # lns_compression = lns_compression_xx_time + lns_compression_yy_time + lns_compression_zz_time + lns_compression_time_2
    # labs_compression = [l.get_label() for l in lns_compression]
    # ax_compression_time.legend(lns_compression, labs_compression, loc=0)
    # ax_compression_time.set_title('Macroscopic stress and strain in compression')
    # # ax_compression_time.set_ylim(ymin=-10,ymax=0)
    # fig_compression_time.tight_layout()
    #
    # fname = fig_dir + 'compression_stress_strain_time'
    # plt.savefig(fname)
    #
    #
    # # ==FIG 5 KINETIC ENERGY IN TENSION=================================================================================
    # fig_KE_tension,ax_KE_tension = plt.subplots()
    # ax_KE_tension.set_ylabel('Kinetic energy [J]')
    # ax_KE_tension.set_xlabel('Time [s]')
    # lns_KE_tension = ax_KE_tension.plot(kinetic_energy_data_tension[: , -1].astype(float),
    #                                     kinetic_energy_data_tension[:,0].astype(float), label=r'KE')
    # ax_KE_tension.set_title('Kinetic energy in tension')
    # plt.legend(loc='best')
    #
    # fname = fig_dir + 'tension_KE'
    # plt.savefig(fname)
    # # ==FIG 6 KINETIC ENERGY IN COMPRESSION=============================================================================
    # fig_KE_compression,ax_KE_compression = plt.subplots()
    # ax_KE_compression.set_ylabel('Kinetic energy [J]')
    # ax_KE_compression.set_xlabel('Time [s]')
    # lns_KE_compression = ax_KE_compression.plot(kinetic_energy_data_compression[: , -1].astype(float), kinetic_energy_data_compression[:,0].astype(float), label=r'KE')
    # ax_KE_compression.set_title('Kinetic energy in compression')
    # plt.legend(loc='best')
    #
    # fname = fig_dir + 'compression_KE'
    # plt.savefig(fname)
    #
    # # ==FIG 7 STIFFNESS AT STRAIN POINTS================================================================================
    if stiffness_at_points_flag == 1:
        strain_points_total_SN_run_1, stiffness_values_total_SN_run_1 = stiffness_func(sxx_tension_SN_run_1,
                                                                                       linear_strain_tension_SN_run_1,
                                                                                       sxx_compression_SN_run_1,
                                                                                       linear_strain_compression_SN_run_1)

        strain_points_total_SN_run_1_El_Pl, stiffness_values_total_SN_run_1_El_Pl = stiffness_func(
            sxx_tension_SN_run_1_El_Pl,
            linear_strain_tension_SN_run_1_El_Pl,
            sxx_compression_SN_run_1_El_Pl,
            linear_strain_compression_SN_run_1_El_Pl)

        strain_points_total_SN_run_1_2E, stiffness_values_total_SN_run_1_2E = stiffness_func(sxx_tension_SN_run_1_2E,
                                                                                             linear_strain_tension_SN_run_1_2E,
                                                                                             sxx_compression_SN_run_1_2E,
                                                                                             linear_strain_compression_SN_run_1_2E)
        strain_points_total_SN_run_1_br_10, stiffness_values_total_SN_run_1_br_10 = stiffness_func(
            sxx_tension_SN_run_1_br_10,
            linear_strain_tension_SN_run_1_br_10,
            sxx_compression_SN_run_1_br_10,
            linear_strain_compression_SN_run_1_br_10)

        strain_points_total_SN_run_1_br_175, stiffness_values_total_SN_run_1_br_175 = stiffness_func(
            sxx_tension_SN_run_1_br_175,
            linear_strain_tension_SN_run_1_br_175,
            sxx_compression_SN_run_1_br_175,
            linear_strain_compression_SN_run_1_br_175)

        # strain_points_total_SN_run_1_EL_particle, stiffness_values_total_SN_run_1_EL_particle = stiffness_func(sxx_tension_SN_run_1_EL_particle,
        #                                                                                linear_strain_tension_SN_run_1_EL_particle,
        #                                                                                sxx_compression_SN_run_1_EL_particle,
        #                                                                                linear_strain_compression_SN_run_1_EL_particle)

        # strain_points_total_SN_run_1_q_0_rmin_3_rmax_10, stiffness_values_total_SN_run_1_q_0_rmin_3_rmax_10 = stiffness_func(sxx_tension_SN_run_1_q_0_rmin_3_rmax_10,
        #                                                                                linear_strain_tension_SN_run_1_q_0_rmin_3_rmax_10,
        #                                                                                sxx_compression_SN_run_1_q_0_rmin_3_rmax_10,
        #                                                                                linear_strain_compression_SN_run_1_q_0_rmin_3_rmax_10)
        #
        # strain_points_total_SN_run_1_q_0_rmin_4_rmax_10, stiffness_values_total_SN_run_1_q_0_rmin_4_rmax_10 = stiffness_func(sxx_tension_SN_run_1_q_0_rmin_4_rmax_10,
        #                                                                                linear_strain_tension_SN_run_1_q_0_rmin_4_rmax_10,
        #                                                                                sxx_compression_SN_run_1_q_0_rmin_4_rmax_10,
        #                                                                                linear_strain_compression_SN_run_1_q_0_rmin_4_rmax_10)

        # strain_points_total_SN_run_1_br_05, stiffness_values_total_SN_run_1_br_05 = stiffness_func(sxx_tension_SN_run_1_br_05,
        #                                                                                linear_strain_tension_SN_run_1_br_05,
        #                                                                                sxx_compression_SN_run_1_br_05,
        #                                                                                linear_strain_compression_SN_run_1_br_05)

        fig_stiff, ax_stiff = plt.subplots()
        ax_stiff.set_ylabel('Unloading stiffness [GPa]')
        ax_stiff.set_xlabel('Strain [%]')
        # ax_stiff.set_title('Unloading stiffness of electrode layer')

        # lns_stiff_EL_particle = ax_stiff.plot(strain_points_total_SN_run_1_EL_particle * 100, stiffness_values_total_SN_run_1_EL_particle * 1E-9,
        #                              linestyle='dashed',
        #                              marker='*', markersize=12, markeredgewidth=3, linewidth=3,
        #                              label=r'EL Particle')

        lns_stiff_2E = ax_stiff.plot(strain_points_total_SN_run_1_2E * 100, stiffness_values_total_SN_run_1_2E * 1E-9,
                                     linestyle='dashed',
                                     marker='+', markersize=12, markeredgewidth=3, linewidth=3,
                                     label=r'$\frac{b_r}{R} = 2.0$')
        lns_stiff_br_175 = ax_stiff.plot(strain_points_total_SN_run_1_br_175 * 100,
                                         stiffness_values_total_SN_run_1_br_175 * 1E-9,
                                         linestyle='dashed',
                                         marker='o', fillstyle='none', markersize=12, markeredgewidth=3, linewidth=3,
                                         label=r'$\frac{b_r}{R} = 1.75$')

        lns_stiff_br_15 = ax_stiff.plot(strain_points_total_SN_run_1 * 100, stiffness_values_total_SN_run_1 * 1E-9,
                                        linestyle='dashed',
                                        marker='x', markersize=12, markeredgewidth=3, linewidth=3,
                                        label=r'$\frac{b_r}{R} = 1.5$')

        lns_stiff_br_10 = ax_stiff.plot(strain_points_total_SN_run_1_br_10 * 100,
                                        stiffness_values_total_SN_run_1_br_10 * 1E-9, linestyle='dashed',
                                        marker='s', fillstyle='none', markersize=12, markeredgewidth=3, linewidth=3,
                                        label=r'$\frac{b_r}{R} = 1.0$')

        # lns_stiff_br_05 = ax_stiff.plot(strain_points_total_SN_run_1_br_05 * 100, stiffness_values_total_SN_run_1_br_05 * 1E-9, linestyle='dashed',
        #                           marker='x', markersize=12, markeredgewidth=3, linewidth=3, label=r'$b_r = 0.5$')
        # lns_label = ax_stiff.plot([], [], ' ', label="Experiments")

        lns_stiff_exp_eps_dot_01 = ax_stiff.plot(exp_strain_points, Modulus_eps_dot_01,
                                                 marker='o', markersize=10, markeredgewidth=3,
                                                 linewidth=0, color='C9', label='$\dot{\Delta} = 0.1 mm/min$')
        lns_stiff_exp_eps_dot_05 = ax_stiff.plot(exp_strain_points, Modulus_eps_dot_05,
                                                 marker='v', markersize=10, markeredgewidth=3,
                                                 linewidth=0, color='C8', label='$\dot{\Delta} = 0.5 mm/min$')
        lns_stiff_exp_eps_dot_10 = ax_stiff.plot(exp_strain_points, Modulus_eps_dot_10,
                                                 marker='s', markersize=10, markeredgewidth=3,
                                                 linewidth=0, color='C7', label='$\dot{\Delta} = 1.0 mm/min$')
        lns_stiff_exp_eps_dot_100 = ax_stiff.plot(exp_strain_points, Modulus_eps_dot_100,
                                                  marker='P', markersize=10, markeredgewidth=3,
                                                  linewidth=0, color='C6', label='$\dot{\Delta} = 10 mm/min$')
        lns_stiff_exp_eps_dot_300 = ax_stiff.plot(exp_strain_points, Modulus_eps_dot_300,
                                                  marker='D', markersize=10, markeredgewidth=3,
                                                  linewidth=0, color='C5', label='$\dot{\Delta} = 30 mm/min$')
        # ===========Legends=================================
        handles_exp = lns_stiff_exp_eps_dot_01 + lns_stiff_exp_eps_dot_05 + lns_stiff_exp_eps_dot_10 + lns_stiff_exp_eps_dot_100 + lns_stiff_exp_eps_dot_300
        handles_exp.reverse()
        labels_exp = [l.get_label() for l in handles_exp]
        first_legend_br = plt.legend(handles_exp, labels_exp, loc='upper right', title='Experimental results')

        ax_stiff.add_artist(first_legend_br)

        handles_br = lns_stiff_2E + lns_stiff_br_175 + lns_stiff_br_15 + lns_stiff_br_10  # +lns_stiff_br_05
        labels_br = [l.get_label() for l in handles_br]
        plt.legend(handles_br, labels_br, loc='upper center', title='Simulations')

        ax_stiff.set_ylim(ymin=0)

        fname = fig_dir + 'stiffness_points_br'
        plt.savefig(fname)

        fig_stiff_El_Pl, ax_stiff_El_Pl = plt.subplots()
        ax_stiff_El_Pl.set_ylabel('Unloading stiffness [GPa]')
        ax_stiff_El_Pl.set_xlabel('Strain [%]')
        # ax_stiff_El_Pl.set_title('Unloading stiffness of electrode layer')

        # lns_stiff_EL_particle = ax_stiff.plot(strain_points_total_SN_run_1_EL_particle * 100, stiffness_values_total_SN_run_1_EL_particle * 1E-9,
        #                              linestyle='dashed',
        #                              marker='*', markersize=12, markeredgewidth=3, linewidth=3,
        #                              label=r'$EL Particle$')

        lns_stiff_El = ax_stiff_El_Pl.plot(strain_points_total_SN_run_1 * 100, stiffness_values_total_SN_run_1 * 1E-9,
                                        linestyle='dashed',
                                        marker='x', markersize=12, markeredgewidth=3, linewidth=3, label=r'$El$')

        lns_stiff_El_Pl = ax_stiff_El_Pl.plot(strain_points_total_SN_run_1_El_Pl * 100,
                                        stiffness_values_total_SN_run_1_El_Pl * 1E-9, linestyle='dashed',
                                        marker='o',fillstyle='none', markersize=12, markeredgewidth=3, linewidth=3, label=r'$El-Pl$')

        # lns_label = ax_stiff.plot([], [], ' ', label="Experiments")

        lns_stiff_exp_eps_dot_01_El_Pl = ax_stiff_El_Pl.plot(exp_strain_points, Modulus_eps_dot_01,
                                                             marker='o', markersize=10, markeredgewidth=3,
                                                             linewidth=0, color='C9',
                                                             label=r'$\dot{\Delta} = 0.1 mm/min$')
        lns_stiff_exp_eps_dot_05_El_Pl = ax_stiff_El_Pl.plot(exp_strain_points, Modulus_eps_dot_05,
                                                             marker='v', markersize=10, markeredgewidth=3,
                                                             linewidth=0, color='C8',
                                                             label='$\dot{\Delta} = 0.5 mm/min$')
        lns_stiff_exp_eps_dot_10_El_Pl = ax_stiff_El_Pl.plot(exp_strain_points, Modulus_eps_dot_10,
                                                             marker='s', markersize=10, markeredgewidth=3,
                                                             linewidth=0, color='C7',
                                                             label='$\dot{\Delta} = 1.0 mm/min$')
        lns_stiff_exp_eps_dot_100_El_Pl = ax_stiff_El_Pl.plot(exp_strain_points, Modulus_eps_dot_100,
                                                              marker='P', markersize=10, markeredgewidth=3,
                                                              linewidth=0, color='C6',
                                                              label='$\dot{\Delta} = 10 mm/min$')
        lns_stiff_exp_eps_dot_300_El_Pl = ax_stiff_El_Pl.plot(exp_strain_points, Modulus_eps_dot_300,
                                                              marker='D', markersize=10, markeredgewidth=3,
                                                              linewidth=0, color='C5',
                                                              label='$\dot{\Delta} = 30 mm/min$')
        # ===========Legends=================================
        handles_exp = lns_stiff_exp_eps_dot_01_El_Pl + lns_stiff_exp_eps_dot_05_El_Pl + lns_stiff_exp_eps_dot_10_El_Pl + lns_stiff_exp_eps_dot_100_El_Pl + lns_stiff_exp_eps_dot_300_El_Pl
        handles_exp.reverse()
        labels_exp = [l.get_label() for l in handles_exp]
        first_legend_El_Pl = plt.legend(handles_exp, labels_exp, loc='upper right', title='Experimental results')

        ax_stiff_El_Pl.add_artist(first_legend_El_Pl)

        handles_El_Pl = lns_stiff_El + lns_stiff_El_Pl
        labels_El_Pl = [l.get_label() for l in handles_El_Pl]
        plt.legend(handles_El_Pl, labels_El_Pl, loc='upper center', title='Simulations')

        ax_stiff_El_Pl.set_ylim(ymin=0)
        fname = fig_dir + 'stiffness_points_El_Pl'
        plt.savefig(fname)

        # ==FIG 8 CONTACTS FOR STRAINS==================================================================================
        """"
        time_vec_compression_SN_run_1, particle_contact_vec_compression_SN_run_1, \
        binder_contact_vec_compression_SN_run_1, binder_particle_contact_vec_compression_SN_run_1 = \
            contact_counter_bertil(simulation_directory_SN_run_1 + '_compression')
        time_vec_compression_SN_run_1_short, linear_strain_compression_SN_run_1_short = time_duplicate_remover(
            time_compression_SN_run_1, linear_strain_compression_SN_run_1)

        time_vec_tension_SN_run_1, particle_contact_vec_tension_SN_run_1, binder_contact_vec_tension_SN_run_1, \
        binder_particle_contact_vec_tension_SN_run_1 = contact_counter_bertil(simulation_directory_SN_run_1
                                                                              + '_tension')
        time_vec_tension_SN_run_1_short, linear_strain_tension_SN_run_1_short = time_duplicate_remover(
            time_tension_SN_run_1, linear_strain_tension_SN_run_1)

        time_vec_compression_SN_run_1_2E, particle_contact_vec_compression_SN_run_1_2E, \
        binder_contact_vec_compression_SN_run_1_2E, binder_particle_contact_vec_compression_SN_run_1_2E = \
            contact_counter_bertil(simulation_directory_SN_run_1_2E + '_compression')
        time_vec_compression_SN_run_1_2E_short, linear_strain_compression_SN_run_1_2E_short = time_duplicate_remover(
            time_compression_SN_run_1_2E, linear_strain_compression_SN_run_1_2E)

        time_vec_tension_SN_run_1_2E, particle_contact_vec_tension_SN_run_1_2E, binder_contact_vec_tension_SN_run_1_2E, \
        binder_particle_contact_vec_tension_SN_run_1_2E = contact_counter_bertil(simulation_directory_SN_run_1_2E
                                                                                 + '_tension')
        time_vec_tension_SN_run_1_2E_short, linear_strain_tension_SN_run_1_2E_short = time_duplicate_remover(
            time_tension_SN_run_1_2E, linear_strain_tension_SN_run_1_2E)

        # time_vec_compression_SN_run_1_br_05, particle_contact_vec_compression_SN_run_1_br_05, \
        # binder_contact_vec_compression_SN_run_1_br_05, binder_particle_contact_vec_compression_SN_run_1_br_05 =\
        #     contact_counter_bertil(simulation_directory_SN_run_1_br_05 + '_compression')
        # time_vec_compression_SN_run_1_br_05_short, linear_strain_compression_SN_run_1_br_05_short= time_duplicate_remover(
        #     time_compression_SN_run_1_br_05,linear_strain_compression_SN_run_1_br_05)
        #
        # time_vec_tension_SN_run_1_br_05, particle_contact_vec_tension_SN_run_1_br_05, binder_contact_vec_tension_SN_run_1_br_05,\
        # binder_particle_contact_vec_tension_SN_run_1_br_05 = contact_counter_bertil(simulation_directory_SN_run_1_br_05
        #                                                                       + '_tension')
        # time_vec_tension_SN_run_1_br_05_short, linear_strain_tension_SN_run_1_br_05_short = time_duplicate_remover(
        #     time_tension_SN_run_1_br_05,linear_strain_tension_SN_run_1_br_05)

        time_vec_compression_SN_run_1_br_10, particle_contact_vec_compression_SN_run_1_br_10, \
        binder_contact_vec_compression_SN_run_1_br_10, binder_particle_contact_vec_compression_SN_run_1_br_10 =\
            contact_counter_bertil(simulation_directory_SN_run_1_br_10 + '_compression')
        time_vec_compression_SN_run_1_br_10_short, linear_strain_compression_SN_run_1_br_10_short= time_duplicate_remover(
            time_compression_SN_run_1_br_10,linear_strain_compression_SN_run_1_br_10)

        time_vec_tension_SN_run_1_br_10, particle_contact_vec_tension_SN_run_1_br_10, binder_contact_vec_tension_SN_run_1_br_10,\
        binder_particle_contact_vec_tension_SN_run_1_br_10 = contact_counter_bertil(simulation_directory_SN_run_1_br_10
                                                                              + '_tension')
        time_vec_tension_SN_run_1_br_10_short, linear_strain_tension_SN_run_1_br_10_short = time_duplicate_remover(
            time_tension_SN_run_1_br_10,linear_strain_tension_SN_run_1_br_10)

        time_vec_compression_SN_run_1_br_175, particle_contact_vec_compression_SN_run_1_br_175, \
        binder_contact_vec_compression_SN_run_1_br_175, binder_particle_contact_vec_compression_SN_run_1_br_175 =\
            contact_counter_bertil(simulation_directory_SN_run_1_br_175 + '_compression')
        time_vec_compression_SN_run_1_br_175_short, linear_strain_compression_SN_run_1_br_175_short= time_duplicate_remover(
            time_compression_SN_run_1_br_175,linear_strain_compression_SN_run_1_br_175)

        time_vec_tension_SN_run_1_br_175, particle_contact_vec_tension_SN_run_1_br_175, binder_contact_vec_tension_SN_run_1_br_175,\
        binder_particle_contact_vec_tension_SN_run_1_br_175 = contact_counter_bertil(simulation_directory_SN_run_1_br_175
                                                                              + '_tension')
        time_vec_tension_SN_run_1_br_175_short, linear_strain_tension_SN_run_1_br_175_short = time_duplicate_remover(
            time_tension_SN_run_1_br_175,linear_strain_tension_SN_run_1_br_175)

        # print(np.where(
        #     np.abs(linear_strain_tension_SN_run_1_short) == np.abs(linear_strain_tension_SN_run_1_short).max())[0][-1])

        # ==PLOTTING DEFINITION======================================================================================
        fig_strain_to_contacts, ax_strain_to_contacts = plt.subplots()

        lns_contacts_tesnion_br_20 = ax_strain_to_contacts.plot(linear_strain_tension_SN_run_1_2E_short[:np.where(
            np.abs(linear_strain_tension_SN_run_1_2E_short) == np.abs(linear_strain_tension_SN_run_1_2E_short).max())[0][
            -1]] * 100, particle_contact_vec_tension_SN_run_1_2E[:np.where(
                                                                        np.abs(
                                                                            linear_strain_tension_SN_run_1_2E_short) == np.abs(
                                                                            linear_strain_tension_SN_run_1_2E_short).max())[
                                                                        0][-1]],color='C0')
        lns_contacts_compression_br_20 = ax_strain_to_contacts.plot(linear_strain_compression_SN_run_1_2E_short[:np.where(
            np.abs(linear_strain_compression_SN_run_1_2E_short) == np.abs(linear_strain_compression_SN_run_1_2E_short).max())[
            0][-1]] * 100,
                                                                particle_contact_vec_compression_SN_run_1_2E[:np.where(
                                                                    np.abs(
                                                                        linear_strain_compression_SN_run_1_2E_short) == np.abs(
                                                                        linear_strain_compression_SN_run_1_2E_short).max())[
                                                                    0][-1]], label=r'$\frac{b_r}{R} = 2.0$',color='C0')

        lns_contacts_compression_br_175 = ax_strain_to_contacts.plot(linear_strain_compression_SN_run_1_br_175_short[:np.where(
            np.abs(linear_strain_compression_SN_run_1_br_175_short) == np.abs(linear_strain_compression_SN_run_1_br_175_short).max())[0][
            -1]] * 100,
                                                                   particle_contact_vec_compression_SN_run_1_br_175[:np.where(
                                                                        np.abs(
                                                                            linear_strain_compression_SN_run_1_br_175_short) == np.abs(
                                                                            linear_strain_compression_SN_run_1_br_175_short).max())[
                                                                        0][-1]],color='C1')
        lns_contacts_tension_br_175 = ax_strain_to_contacts.plot(linear_strain_tension_SN_run_1_br_175_short[:np.where(
            np.abs(linear_strain_tension_SN_run_1_br_175_short) == np.abs(linear_strain_tension_SN_run_1_br_175_short).max())[
            0][-1]] * 100,
                                                                particle_contact_vec_tension_SN_run_1_br_175[:np.where(
                                                                    np.abs(
                                                                        linear_strain_tension_SN_run_1_br_175_short) == np.abs(
                                                                        linear_strain_tension_SN_run_1_br_175_short).max())[
                                                                    0][-1]], label=r'$\frac{b_r}{R} = 1.75$',color='C1')


        lns_contacts_compression_br_15 = ax_strain_to_contacts.plot(linear_strain_tension_SN_run_1_short[:np.where(
            np.abs(linear_strain_tension_SN_run_1_short) == np.abs(linear_strain_tension_SN_run_1_short).max())[0][
            -1]] * 100,
                                                                    particle_contact_vec_tension_SN_run_1[:np.where(
                                                                        np.abs(
                                                                            linear_strain_tension_SN_run_1_short) == np.abs(
                                                                            linear_strain_tension_SN_run_1_short).max())[
                                                                        0][-1]],color='C2')
        lns_contacts_tension_br_15 = ax_strain_to_contacts.plot(linear_strain_compression_SN_run_1_short[:np.where(
            np.abs(linear_strain_compression_SN_run_1_short) == np.abs(linear_strain_compression_SN_run_1_short).max())[
            0][-1]] * 100,
                                                                particle_contact_vec_compression_SN_run_1[:np.where(
                                                                    np.abs(
                                                                        linear_strain_compression_SN_run_1_short) == np.abs(
                                                                        linear_strain_compression_SN_run_1_short).max())[
                                                                    0][-1]], label=r'$\frac{b_r}{R} = 1.5$',color='C2')

        lns_contacts_compression_br_10 = ax_strain_to_contacts.plot(linear_strain_compression_SN_run_1_br_10_short[:np.where(
            np.abs(linear_strain_compression_SN_run_1_br_10_short) == np.abs(linear_strain_compression_SN_run_1_br_10_short).max())[0][
            -1]] * 100,
                                                                    particle_contact_vec_compression_SN_run_1_br_10[:np.where(
                                                                        np.abs(
                                                                            linear_strain_compression_SN_run_1_br_10_short) == np.abs(
                                                                            linear_strain_compression_SN_run_1_br_10_short).max())[
                                                                        0][-1]],color='C3')
        lns_contacts_tension_br_10 = ax_strain_to_contacts.plot(linear_strain_tension_SN_run_1_br_10_short[:np.where(
            np.abs(linear_strain_tension_SN_run_1_br_10_short) == np.abs(linear_strain_tension_SN_run_1_br_10_short).max())[
            0][-1]] * 100,
                                                                particle_contact_vec_tension_SN_run_1_br_10[:np.where(
                                                                    np.abs(
                                                                        linear_strain_tension_SN_run_1_br_10_short) == np.abs(
                                                                        linear_strain_tension_SN_run_1_br_10_short).max())[
                                                                    0][-1]], label=r'$\frac{b_r}{R} = 1.0$',color='C3')


        # lns_contacts_compression_br_05 = ax_strain_to_contacts.plot(linear_strain_compression_SN_run_1_br_05_short[:np.where(
        #     np.abs(linear_strain_compression_SN_run_1_br_05_short) == np.abs(linear_strain_compression_SN_run_1_br_05_short).max())[0][
        #     -1]] * 100,
        #                                                             particle_contact_vec_compression_SN_run_1_br_05[:np.where(
        #                                                                 np.abs(
        #                                                                     linear_strain_compression_SN_run_1_br_05_short) == np.abs(
        #                                                                     linear_strain_compression_SN_run_1_br_05_short).max())[
        #                                                                 0][-1]],color='C4')
        # lns_contacts_tension_br_05 = ax_strain_to_contacts.plot(linear_strain_tension_SN_run_1_br_05_short[:np.where(
        #     np.abs(linear_strain_tension_SN_run_1_br_05_short) == np.abs(linear_strain_tension_SN_run_1_br_05_short).max())[
        #     0][-1]] * 100,
        #                                                         particle_contact_vec_tension_SN_run_1_br_05[:np.where(
        #                                                             np.abs(
        #                                                                 linear_strain_tension_SN_run_1_br_05_short) == np.abs(
        #                                                                 linear_strain_tension_SN_run_1_br_05_short).max())[
        #                                                             0][-1]], label=r'$\frac{b_r}{R} = 0.5$',color='C4')


        ax_strain_to_contacts.set_ylim(ymin=0)
        ax_strain_to_contacts.set_xlim(xmin=-2.2, xmax=2.2)
        ax_strain_to_contacts.set_ylabel('Particle contacts [-]')
        ax_strain_to_contacts.set_xlabel('Strain [%]')
        ax_strain_to_contacts.legend(loc='best')
        fname = fig_dir + 'contacts_to_strain'
        plt.savefig(fname)
        """""
        # ==============================================================================================================

        # # =FIGURE COMPARING SIZE DISTRIBUTIONS========================================================================
        # fig_stiff_size_distr, ax_stiff_size_distr = plt.subplots()
        # ax_stiff_size_distr.set_ylabel('Unloading stiffness [GPa]')
        # ax_stiff_size_distr.set_xlabel('Strain [%]')
        # # ax_stiff_El_Pl.set_title('Unloading stiffness of electrode layer')
        #
        # lns_stiff_reference = ax_stiff_size_distr.plot(strain_points_total_SN_run_1 * 100, stiffness_values_total_SN_run_1 * 1E-9,
        #                                    linestyle='dashed',
        #                                    marker='x', markersize=12, markeredgewidth=3, linewidth=3, label=r'$El$')
        #
        # lns_stiff_q_0_rmin_3 = ax_stiff_size_distr.plot(strain_points_total_SN_run_1_q_0_rmin_3_rmax_10 * 100,
        #                                       stiffness_values_total_SN_run_1_q_0_rmin_3_rmax_10 * 1E-9, linestyle='dashed',
        #                                       marker='o', fillstyle='none', markersize=12, markeredgewidth=3, linewidth=3,
        #                                       label=r'$r_min_3$')
        # lns_stiff_q_0_rmin_4 = ax_stiff_size_distr.plot(strain_points_total_SN_run_1_q_0_rmin_4_rmax_10 * 100,
        #                                       stiffness_values_total_SN_run_1_q_0_rmin_4_rmax_10 * 1E-9, linestyle='dashed',
        #                                       marker='*', fillstyle='none', markersize=12, markeredgewidth=3, linewidth=3,
        #                                       label=r'$r_min_4$')
        #
        # # lns_label = ax_stiff.plot([], [], ' ', label="Experiments")
        #
        # lns_stiff_exp_eps_dot_01_El_Pl = ax_stiff_size_distr.plot(exp_strain_points, Modulus_eps_dot_01,
        #                                                      marker='o', markersize=10, markeredgewidth=3,
        #                                                      linewidth=0, color='C9', label=r'$\dot{\Delta} = 0.1 mm/min$')
        # lns_stiff_exp_eps_dot_05_El_Pl = ax_stiff_size_distr.plot(exp_strain_points, Modulus_eps_dot_05,
        #                                                      marker='v', markersize=10, markeredgewidth=3,
        #                                                      linewidth=0, color='C8', label='$\dot{\Delta} = 0.5 mm/min$')
        # lns_stiff_exp_eps_dot_10_El_Pl = ax_stiff_size_distr.plot(exp_strain_points, Modulus_eps_dot_10,
        #                                                      marker='s', markersize=10, markeredgewidth=3,
        #                                                      linewidth=0, color='C7', label='$\dot{\Delta} = 1.0 mm/min$')
        # lns_stiff_exp_eps_dot_100_El_Pl = ax_stiff_size_distr.plot(exp_strain_points, Modulus_eps_dot_100,
        #                                                       marker='P', markersize=10, markeredgewidth=3,
        #                                                       linewidth=0, color='C6', label='$\dot{\Delta} = 10 mm/min$')
        # lns_stiff_exp_eps_dot_300_El_Pl = ax_stiff_size_distr.plot(exp_strain_points, Modulus_eps_dot_300,
        #                                                       marker='D', markersize=10, markeredgewidth=3,
        #                                                       linewidth=0, color='C5', label='$\dot{\Delta} = 30 mm/min$')
        # # ===========Legends=================================
        # handles_exp = lns_stiff_exp_eps_dot_01_El_Pl + lns_stiff_exp_eps_dot_05_El_Pl + lns_stiff_exp_eps_dot_10_El_Pl + lns_stiff_exp_eps_dot_100_El_Pl + lns_stiff_exp_eps_dot_300_El_Pl
        # labels_exp = [l.get_label() for l in handles_exp]
        # first_legend_size_distr = plt.legend(handles_exp, labels_exp, loc='upper right', title='Experimental results')
        #
        # ax_stiff_size_distr.add_artist(first_legend_size_distr)
        #
        # handles_size_distr = lns_stiff_reference + lns_stiff_q_0_rmin_3 + lns_stiff_q_0_rmin_4
        # labels_size_distr = [l.get_label() for l in handles_size_distr]
        # plt.legend(handles_size_distr, labels_size_distr, loc='upper center', title='Simulations')
        #
        # ax_stiff_size_distr.set_ylim(ymin=0)
        # fname = fig_dir + 'stiffness_points_size_distr'
        # plt.savefig(fname)

    # =HISTOGRAM OF PARTICLE DISTRIBUTIONS==============================================================================
    # reference_particle_data = one_file_reader(
    #     simulation_directory_SN_run_1 + '_compression/particles/particles_' + str(time_tension_SN_run_1[-1]) + '.dou')
    # reference_particle_radius = np.array([])
    # for x in reference_particle_data:
    #     reference_particle_radius = np.append(reference_particle_radius, (float(x.split(", ")[7])))
    #
    #
    # q_0_rmin_3_particle_data = one_file_reader(
    #     simulation_directory_SN_run_1_q_0_rmin_3_rmax_10 + '_compression/particles/particles_' + str(
    #         time_tension_SN_run_1_q_0_rmin_3_rmax_10[-1]) + '.dou')
    # q_0_rmin_3_particle_radius = np.array([])
    # for x in q_0_rmin_3_particle_data:
    #     q_0_rmin_3_particle_radius = np.append(q_0_rmin_3_particle_radius, (float(x.split(", ")[7])))
    #
    #
    # q_0_rmin_4_particle_data = one_file_reader(
    #     simulation_directory_SN_run_1_q_0_rmin_4_rmax_10 + '_compression/particles/particles_' + str(
    #         time_tension_SN_run_1_q_0_rmin_4_rmax_10[-1]) + '.dou')
    # q_0_rmin_4_particle_radius = np.array([])
    # for x in q_0_rmin_4_particle_data:
    #     q_0_rmin_4_particle_radius = np.append(q_0_rmin_4_particle_radius, (float(x.split(", ")[7])))
    #
    #
    # Ebner_particle_data_df = pd.read_csv(
    # "G:/My Drive/Skola/KTH/PhD/Litteratur/X-Ray Tomography of Porous, Transition Metal Oxide Based Lithium Ion Battery Electrodes/NMC_96wt_0bar.dat")
    # Ebner_particle_volume=Ebner_particle_data_df['volume'].to_numpy()
    # Ebner_particle_radius = (Ebner_particle_volume*3/(4*3.14))**(1/3)
    # # Ebner_particle_radius = Ebner_particle_radius[(Ebner_particle_radius > 3) & (Ebner_particle_radius < 10)]
    #
    # # s_q_0[(s_q_0 > D_min) & (s_q_0 < D_max)]
    # # q_0_rmin_4_particle_radius = np.array([])
    # # for x in q_0_rmin_4_particle_data:
    # #     q_0_rmin_4_particle_radius = np.append(q_0_rmin_4_particle_radius, (float(x.split(", ")[7])))
    #
    # fig_size_hist,ax_size_hist = plt.subplots(1,4)
    # reference_lns = ax_size_hist[0].hist(reference_particle_radius*100,bins=50,density=True)#bins=50, label='ref')
    # q_0_rmin_3_lns = ax_size_hist[1].hist(q_0_rmin_3_particle_radius*100,bins=50,density=True)#bins=50, label='r_min=3µm')
    # q_0_rmin_4_lns = ax_size_hist[2].hist(q_0_rmin_4_particle_radius*100,bins=50,density=True)#bins=50, label='r_min=4µm')
    # Ebner_lns = ax_size_hist[3].hist(Ebner_particle_radius,bins=50,density=True)#bins=50)
    # ax_size_hist[0].set_title('ref')
    # ax_size_hist[1].set_title('r_min=3µm')
    # ax_size_hist[2].set_title('r_min=4µm')
    # ax_size_hist[3].set_title('Ebner')
    # for x in ax_size_hist:
    #     x.set_xlabel('Particle radius [µm]')
    #     x.set_ylabel('Number of particles')
    #     # x.set_ylim(ymax=400)
    #     x.set_xlim(xmin=3,xmax=12)
    #     # ax_size_hist.legend(loc='best')
    # ==SHOW PLOT=======================================================================================================
    plt.show()
