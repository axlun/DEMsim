from Bertil_functions.Bertil_functions import *

import numpy as np
import matplotlib
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


def elastic_springback_func(calendering_surface_pressure,calendering_surface_position):
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
    # print(loading_height)
    #
    # print(min_height)
    #
    # print(unload_height)
    #
    print(specific_elastic_springback)
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
    fig_dir = 'C:/temp/figures/Bertil_calendering_pressure_multiple_simulation/'
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

    calendering_time_SN_run_1, calendering_surface_pressure_SN_run_1, bottom_surface_pressure_SN_run_1, calendering_surface_position_SN_run_1, sxx_SN_run_1, syy_SN_run_1, szz_SN_run_1, tau_xy_SN_run_1, tau_xz_SN_run_1, tau_yx_SN_run_1, tau_yz_SN_run_1, tau_zx_SN_run_1, tau_zy_SN_run_1,kinetic_energy_data_SN_run_1 = calendering_plot_processing(simulation_directory_SN_run_1)

    calendering_time_SN_run_1_spread_1, calendering_surface_pressure_SN_run_1_spread_1, bottom_surface_pressure_SN_run_1_spread_1, calendering_surface_position_SN_run_1_spread_1, sxx_SN_run_1_spread_1, syy_SN_run_1_spread_1, szz_SN_run_1_spread_1, tau_xy_SN_run_1_spread_1, tau_xz_SN_run_1_spread_1, tau_yx_SN_run_1_spread_1, tau_yz_SN_run_1_spread_1, tau_zx_SN_run_1_spread_1, tau_zy_SN_run_1_spread_1, kinetic_energy_data_SN_run_1_spread_1 = calendering_plot_processing(
        simulation_directory_SN_run_1_spread_1)

    calendering_time_SN_run_1_spread_2, calendering_surface_pressure_SN_run_1_spread_2, bottom_surface_pressure_SN_run_1_spread_2, calendering_surface_position_SN_run_1_spread_2, sxx_SN_run_1_spread_2, syy_SN_run_1_spread_2, szz_SN_run_1_spread_2, tau_xy_SN_run_1_spread_2, tau_xz_SN_run_1_spread_2, tau_yx_SN_run_1_spread_2, tau_yz_SN_run_1_spread_2, tau_zx_SN_run_1_spread_2, tau_zy_SN_run_1_spread_2, kinetic_energy_data_SN_run_1_spread_2 = calendering_plot_processing(
        simulation_directory_SN_run_1_spread_2)

    calendering_time_SN_run_1_spread_3, calendering_surface_pressure_SN_run_1_spread_3, bottom_surface_pressure_SN_run_1_spread_3, calendering_surface_position_SN_run_1_spread_3, sxx_SN_run_1_spread_3, syy_SN_run_1_spread_3, szz_SN_run_1_spread_3, tau_xy_SN_run_1_spread_3, tau_xz_SN_run_1_spread_3, tau_yx_SN_run_1_spread_3, tau_yz_SN_run_1_spread_3, tau_zx_SN_run_1_spread_3, tau_zy_SN_run_1_spread_3, kinetic_energy_data_SN_run_1_spread_3 = calendering_plot_processing(
        simulation_directory_SN_run_1_spread_3)

    calendering_time_SN_run_1_spread_4, calendering_surface_pressure_SN_run_1_spread_4, bottom_surface_pressure_SN_run_1_spread_4, calendering_surface_position_SN_run_1_spread_4, sxx_SN_run_1_spread_4, syy_SN_run_1_spread_4, szz_SN_run_1_spread_4, tau_xy_SN_run_1_spread_4, tau_xz_SN_run_1_spread_4, tau_yx_SN_run_1_spread_4, tau_yz_SN_run_1_spread_4, tau_zx_SN_run_1_spread_4, tau_zy_SN_run_1_spread_4, kinetic_energy_data_SN_run_1_spread_4 = calendering_plot_processing(
        simulation_directory_SN_run_1_spread_4)

    calendering_time_SN_run_1_El_Pl, calendering_surface_pressure_SN_run_1_El_Pl, bottom_surface_pressure_SN_run_1_El_Pl, calendering_surface_position_SN_run_1_El_Pl, sxx_SN_run_1_El_Pl, syy_SN_run_1_El_Pl, szz_SN_run_1_El_Pl, tau_xy_SN_run_1_El_Pl, tau_xz_SN_run_1_El_Pl, tau_yx_SN_run_1_El_Pl, tau_yz_SN_run_1_El_Pl, tau_zx_SN_run_1_El_Pl, tau_zy_SN_run_1_El_Pl, kinetic_energy_data_SN_run_1_El_Pl = calendering_plot_processing(simulation_directory_SN_run_1_El_Pl)
    calendering_time_SN_run_1_br_175, calendering_surface_pressure_SN_run_1_br_175, bottom_surface_pressure_SN_run_1_br_175, calendering_surface_position_SN_run_1_br_175, sxx_SN_run_1_br_175, syy_SN_run_1_br_175, szz_SN_run_1_br_175, tau_xy_SN_run_1_br_175, tau_xz_SN_run_1_br_175, tau_yx_SN_run_1_br_175, tau_yz_SN_run_1_br_175, tau_zx_SN_run_1_br_175, tau_zy_SN_run_1_br_175, kinetic_energy_data_SN_run_1_br_175 = calendering_plot_processing(
        simulation_directory_SN_run_1_br_175)
    calendering_time_SN_run_1_br_10, calendering_surface_pressure_SN_run_1_br_10, bottom_surface_pressure_SN_run_1_br_10, calendering_surface_position_SN_run_1_br_10, sxx_SN_run_1_br_10, syy_SN_run_1_br_10, szz_SN_run_1_br_10, tau_xy_SN_run_1_br_10, tau_xz_SN_run_1_br_10, tau_yx_SN_run_1_br_10, tau_yz_SN_run_1_br_10, tau_zx_SN_run_1_br_10, tau_zy_SN_run_1_br_10, kinetic_energy_data_SN_run_1_br_10 = calendering_plot_processing(simulation_directory_SN_run_1_br_10)
    calendering_time_SN_run_1_br_05, calendering_surface_pressure_SN_run_1_br_05, bottom_surface_pressure_SN_run_1_br_05, calendering_surface_position_SN_run_1_br_05, sxx_SN_run_1_br_05, syy_SN_run_1_br_05, szz_SN_run_1_br_05, tau_xy_SN_run_1_br_05, tau_xz_SN_run_1_br_05, tau_yx_SN_run_1_br_05, tau_yz_SN_run_1_br_05, tau_zx_SN_run_1_br_05, tau_zy_SN_run_1_br_05, kinetic_energy_data_SN_run_1_br_05 = calendering_plot_processing(simulation_directory_SN_run_1_br_05)
    calendering_time_SN_run_1_2E, calendering_surface_pressure_SN_run_1_2E, bottom_surface_pressure_SN_run_1_2E, calendering_surface_position_SN_run_1_2E, sxx_SN_run_1_2E, syy_SN_run_1_2E, szz_SN_run_1_2E, tau_xy_SN_run_1_2E, tau_xz_SN_run_1_2E, tau_yx_SN_run_1_2E, tau_yz_SN_run_1_2E, tau_zx_SN_run_1_2E, tau_zy_SN_run_1_2E, kinetic_energy_data_SN_run_1_2E = calendering_plot_processing(
        simulation_directory_SN_run_1_2E)
    print('Elastic springback with different binder radii')
    loading_height_SN_run_1_2E, unload_height_SN_run_1_2E, specific_elastic_springback_SN_run_1_2E = \
        elastic_springback_func(
         calendering_surface_pressure_SN_run_1_2E,
         calendering_surface_position_SN_run_1_2E)

    loading_height_SN_run_1_br_175, unload_height_SN_run_1_br_175, specific_elastic_springback_SN_run_1_br_175 = \
        elastic_springback_func(
         calendering_surface_pressure_SN_run_1_br_175,
         calendering_surface_position_SN_run_1_br_175)

    loading_height_SN_run_1, unload_height_SN_run_1, specific_elastic_springback_SN_run_1 = \
        elastic_springback_func(
         calendering_surface_pressure_SN_run_1,
         calendering_surface_position_SN_run_1)

    loading_height_SN_run_1_br_10, unload_height_SN_run_1_br_10, specific_elastic_springback_SN_run_1_br_10 = \
        elastic_springback_func(
         calendering_surface_pressure_SN_run_1_br_10,
         calendering_surface_position_SN_run_1_br_10)


    print('El, El-Pl elastic springback')
    elastic_springback_func(calendering_surface_pressure_SN_run_1, calendering_surface_position_SN_run_1)

    loading_height_SN_run_1_El_Pl, unload_height_SN_run_1_El_Pl, specific_elastic_springback_SN_run_1_El_Pl = \
        elastic_springback_func(
         calendering_surface_pressure_SN_run_1_El_Pl,
         calendering_surface_position_SN_run_1_El_Pl)






    # # ===FIG 1 CALENDERING SURFACE PRESSURE SURFACE POSITION TIME=======================================================
    # fig_calendering_surface_pressure, ax_calendering_surface_pressure = plt.subplots()
    # ax_calendering_surface_pressure.set_ylabel("Calendering surface pressure [MPa]")
    # ax_calendering_surface_pressure.set_xlabel("Time [s]")
    # lns_calendering_surface_pressure = ax_calendering_surface_pressure.plot(calendering_time_SN_run_1,
    #                                                                         calendering_surface_pressure_SN_run_1 * 1E-6,
    #                                                                         label='Pressure')
    # ax_calendering_surface_position = ax_calendering_surface_pressure.twinx()
    # ax_calendering_surface_position.set_ylabel("Calendering surface height [m]")
    # lns_calendering_surface_position = ax_calendering_surface_position.plot(calendering_time_SN_run_1,
    #                                                                         calendering_surface_position_SN_run_1, 'c',
    #                                                                         label='Position')
    # # ax_calendering_surface_pressure.set_title('Calendering surface pressure')
    #
    # lns = lns_calendering_surface_pressure + lns_calendering_surface_position
    # labs = [l.get_label() for l in lns]
    # ax_calendering_surface_pressure.legend(lns, labs, loc='best')
    # fig_calendering_surface_pressure.tight_layout()
    #
    # fname = fig_dir + 'calendering_surface_pressure_surface_position_time'
    # plt.savefig(fname)
    #
    # # ===FIG 2 CALENDERING SURFACE PRESSURE POSITION====================================================================
    # calendering_initiate_index, calendering_break_index = calendering_break_index_func(calendering_surface_position_SN_run_1)
    # fig_calendering_process, ax_calendering_process = plt.subplots()
    # ax_calendering_process.plot(calendering_surface_position_SN_run_1[calendering_initiate_index:(calendering_break_index - 1)],
    #                             calendering_surface_pressure_SN_run_1[
    #                             calendering_initiate_index:(calendering_break_index - 1)] * 1e-6)
    # ax_calendering_process.set_ylabel("Calendering surface pressure [MPa]")
    # ax_calendering_process.set_xlabel("Calendering surface height [m]")
    # # ax_calendering_process.set_title('Calendering surface pressure')
    # fig_calendering_process.tight_layout()
    #
    # fname = fig_dir + 'calendering_surface_pressure_surface_position'
    # plt.savefig(fname)

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

    ax_full_calendering_sim.set_xlim(xmin=104.8,xmax=120)
    ax_full_calendering_sim.xaxis.set_major_locator(MultipleLocator(5))

    ax_full_calendering_sim.set_ylim(ymin=0)
    ax_full_calendering_sim.set_ylabel("Calendering surface pressure [MPa]")
    ax_full_calendering_sim.set_xlabel("Calendering surface height [µm]")
    #    ax_full_calendering_sim.set_title('calendering surface pressure')
    fig_full_calendering_sim.tight_layout()
    ax_full_calendering_sim.legend(loc="best")

    fname = fig_dir + 'calendering_surface_pressure_surface_position_full_statistical_spread'
    plt.savefig(fname)

    # ===FIG 3 EL/ EL-PL CALENDERING SURFACE PRESSURE SURFACE POSITION FULL=======================================================
    fig_full_calendering_sim, ax_full_calendering_sim = plt.subplots()
    ax_full_calendering_sim.plot(calendering_surface_position_SN_run_1[:] * 1E2, calendering_surface_pressure_SN_run_1[:] * 1E-6, label=r'$El$')
    ax_full_calendering_sim.plot(calendering_surface_position_SN_run_1_El_Pl[:] * 1E2,
                                 calendering_surface_pressure_SN_run_1_El_Pl[:] * 1E-6, label=r'$El-Pl$')
    ax_full_calendering_sim.set_xlim(xmin=104.8,xmax=120)
    ax_full_calendering_sim.xaxis.set_major_locator(MultipleLocator(5))

    ax_full_calendering_sim.set_ylim(ymin=0)
    ax_full_calendering_sim.set_ylabel("Calendering surface pressure [MPa]")
    ax_full_calendering_sim.set_xlabel("Calendering surface height [µm]")
    #    ax_full_calendering_sim.set_title('calendering surface pressure')
    fig_full_calendering_sim.tight_layout()
    ax_full_calendering_sim.legend(loc="best")

    fname = fig_dir + 'calendering_surface_pressure_surface_position_full_El_Pl'
    plt.savefig(fname)

    # ===FIG  CALENDERING SURFACE PRESSURE SURFACE POSITION FULL=======================================================
    fig_full_calendering_sim, ax_full_calendering_sim = plt.subplots()
    ax_full_calendering_sim.plot(calendering_surface_position_SN_run_1_2E[:] * 1E2,
                                 calendering_surface_pressure_SN_run_1_2E[:] * 1E-6, label=r'$\frac{b_r}{R} = 2.0$')
    ax_full_calendering_sim.plot(calendering_surface_position_SN_run_1_br_175[:] * 1E2,
                                 calendering_surface_pressure_SN_run_1_br_175[:] * 1E-6, label=r'$\frac{b_r}{R} = 1.75$')
    ax_full_calendering_sim.plot(calendering_surface_position_SN_run_1[:] * 1E2,
                                 calendering_surface_pressure_SN_run_1[:] * 1E-6, label=r'$\frac{b_r}{R} = 1.5$')
    ax_full_calendering_sim.plot(calendering_surface_position_SN_run_1_br_10[:] * 1E2,
                                 calendering_surface_pressure_SN_run_1_br_10[:] * 1E-6, label=r'$\frac{b_r}{R} = 1.0$')
    # ax_full_calendering_sim.plot(calendering_surface_position_SN_run_1_br_05[:],
    #                              calendering_surface_pressure_SN_run_1_br_05[:] * 1E-6, label=r'$b_r = 0.5$')
    ax_full_calendering_sim.set_xlim(xmin=104.8,xmax=120)
    ax_full_calendering_sim.xaxis.set_major_locator(MultipleLocator(5))

    ax_full_calendering_sim.set_ylim(ymin=0)
    ax_full_calendering_sim.set_ylabel("Calendering surface pressure [MPa]")
    ax_full_calendering_sim.set_xlabel("Calendering surface height [µm]")
    #    ax_full_calendering_sim.set_title('calendering surface pressure')
    fig_full_calendering_sim.tight_layout()
    ax_full_calendering_sim.legend(loc="best")

    fname = fig_dir + 'calendering_surface_pressure_surface_position_full_br'
    plt.savefig(fname)

    # # ===FIG 4 STRESS ZZ SURFACE POSITION===============================================================================
    #
    # fig_sigma_zz_calendering_surface_position, ax_sigma_zz_calendering_surface_position = plt.subplots()
    # ax_sigma_zz_calendering_surface_position.plot(calendering_surface_position_SN_run_1[:], szz_SN_run_1)
    # ax_sigma_zz_calendering_surface_position.set_ylabel("Stress in z [Pa]")
    # ax_sigma_zz_calendering_surface_position.set_xlabel("Calendering surface height [m]")
    # # ax_sigma_zz_calendering_surface_position.set_title("Stress in Z to calendering surface position")
    # fig_sigma_zz_calendering_surface_position.tight_layout()
    #
    # # ===FIG 5 STRESS ZZ SURFACE PRESSURE SURFACE POSITION TIME=========================================================
    # fig_sigma_zz_time, ax_sigma_zz_time = plt.subplots()
    # lns_sigma_zz_time = ax_sigma_zz_time.plot(calendering_time_SN_run_1, -szz_SN_run_1 * 1e-6, label=r'$\sigma_{zz}$')
    # lsn_cal_surf_pre_time = ax_sigma_zz_time.plot(calendering_time_SN_run_1, calendering_surface_pressure_SN_run_1 * 1e-6,
    #                                               label=r'Surface pressure')
    # ax_sigma_zz_time.set_ylabel("Stress in z [MPa]")
    # ax_sigma_zz_time.set_xlabel("Time [s]")
    # ax_calendering_surface_position_time = ax_sigma_zz_time.twinx()
    # ax_calendering_surface_position_time.set_ylabel("Calendering surface height [m]")
    # lns_calendering_surface_position_time = ax_calendering_surface_position_time.plot(calendering_time_SN_run_1,
    #                                                                                   calendering_surface_position_SN_run_1, 'c',
    #                                                                                   label='Surface position')
    # lns2 = lns_sigma_zz_time + lsn_cal_surf_pre_time + lns_calendering_surface_position_time
    # labs2 = [l.get_label() for l in lns2]
    # ax_sigma_zz_time.legend(lns2, labs2, loc="best")
    # fig_sigma_zz_time.tight_layout()
    #
    # # ===FIG 6 STRESS ZZ SURFACE PRESSURE TIME==========================================================================
    # fig_calendering_surface_pressure_sigma_zz, ax_sigma_zz_time2 = plt.subplots()
    # lns5 = ax_sigma_zz_time2.plot(calendering_time_SN_run_1, -szz_SN_run_1, label=r'$\sigma_{zz}$')
    # ax_sigma_zz_time2.set_ylabel("Stress [Pa]")
    # ax_sigma_zz_time2.set_xlabel("Time [s]")
    # lns6 = ax_sigma_zz_time2.plot(calendering_time_SN_run_1, calendering_surface_pressure_SN_run_1[:], label='Surface pressure')
    # #  ax9.set_title('calendering surface pressure')
    # lns3 = lns5 + lns6
    # labs3 = [l.get_label() for l in lns3]
    # ax_sigma_zz_time2.legend(lns3, labs3, loc="best")
    # fig_calendering_surface_pressure_sigma_zz.tight_layout()
    #
    # # ===FIG 7 ALL SURFACE PRESSURES STRESS ZZ SURFACE POSITION TIME====================================================
    # fig_calendering_surface_bottom_surface, ax_calendering_surface_bottom_surface = plt.subplots()
    # lns_calendering_surface_pressure_2 = ax_calendering_surface_bottom_surface.plot(calendering_time_SN_run_1,
    #                                                                                 calendering_surface_pressure_SN_run_1 * 1e-6,
    #                                                                                 'g', label=r'Surface pressure')
    # lns_bottom_surface_pressure = ax_calendering_surface_bottom_surface.plot(calendering_time_SN_run_1,
    #                                                                          bottom_surface_pressure_SN_run_1 * 1e-6,
    #                                                                          label=r'Bottom pressure')
    # lns_macroscopic_stress = ax_calendering_surface_bottom_surface.plot(calendering_time_SN_run_1, -szz_SN_run_1 * 1e-6,
    #                                                                     label=r'$\sigma_{zz}$')
    #
    # ax_calendering_surface_position_2 = ax_calendering_surface_bottom_surface.twinx()
    # ax_calendering_surface_position_2.set_ylabel("Calendering surface height [m]")
    # #    ax_calendering_surface_position_2.set_ylim([0,200])
    #
    # lns_calendering_surface_position_2 = ax_calendering_surface_position_2.plot(calendering_time_SN_run_1,
    #                                                                             calendering_surface_position_SN_run_1, 'c',
    #                                                                             label='Surface position')
    # ax_calendering_surface_bottom_surface.set_ylabel("Stress in z [MPa]")
    # ax_calendering_surface_bottom_surface.set_xlabel("Time [s]")
    # # ax_calendering_surface_bottom_surface.set_title("Pressure on calendering surface and bottom surface")
    # lns_sum = lns_calendering_surface_pressure_2 + lns_bottom_surface_pressure + lns_macroscopic_stress + lns_calendering_surface_position_2
    # labs4 = [l.get_label() for l in lns_sum]
    # # ax_calendering_surface_bottom_surface.set_ylim(ymin=-.02,ymax=.12)
    # ax_calendering_surface_bottom_surface.legend(lns_sum, labs4, loc="best")
    # fig_calendering_surface_bottom_surface.tight_layout()
    #
    # fname = fig_dir + 'all_surface_pressure_stess_ZZ_surface_position_time'
    # plt.savefig(fname)
    #
    # # ===FIG 8 KINETIC  ENERGY===============================================================================================
    # fig_KE, ax_KE = plt.subplots()
    # lns5 = ax_KE.plot(kinetic_energy_data_SN_run_1[:, -1], kinetic_energy_data_SN_run_1[:, 0], label=r'KE')
    # ax_KE.set_ylabel("Kinetic energy [J]")
    # ax_KE.set_xlabel("Time [s]")
    # # ax_KE.set_title("Kinetic energy of particles")
    # ax_KE.legend(loc="best")
    # fig_KE.tight_layout()
    #
    # fname = fig_dir + 'KE'
    # plt.savefig(fname)
    #
    # # ===FIG 9 STRESS XX YY ZZ TIME===============================================================================================
    # fig_force_fabric_stress_time, ax_force_fabric_stress_time = plt.subplots()
    # lns_sxx = ax_force_fabric_stress_time.plot(calendering_time_SN_run_1, -sxx_SN_run_1 * 1e-6, label=r'$\sigma_{xx}$')
    # lns_syy = ax_force_fabric_stress_time.plot(calendering_time_SN_run_1, -syy_SN_run_1 * 1e-6, label=r'$\sigma_{yy}$')
    # lns_szz = ax_force_fabric_stress_time.plot(calendering_time_SN_run_1, -szz_SN_run_1 * 1e-6, label=r'$\sigma_{zz}$')
    # lns_tau_xy = ax_force_fabric_stress_time.plot(calendering_time_SN_run_1, -tau_xy_SN_run_1 * 1e-6, label=r'$\tau_{xy}$')
    # lns_tau_xz = ax_force_fabric_stress_time.plot(calendering_time_SN_run_1, -tau_xz_SN_run_1 * 1e-6, label=r'$\tau_{xz}$')
    # lns_tau_yx = ax_force_fabric_stress_time.plot(calendering_time_SN_run_1, -tau_yx_SN_run_1 * 1e-6, label=r'$\tau_{yx}$')
    # lns_tau_yz = ax_force_fabric_stress_time.plot(calendering_time_SN_run_1, -tau_yz_SN_run_1 * 1e-6, label=r'$\tau_{yz}$')
    # lns_tau_zx = ax_force_fabric_stress_time.plot(calendering_time_SN_run_1, -tau_zx_SN_run_1 * 1e-6, label=r'$\tau_{zx}$')
    # lns_tau_zy = ax_force_fabric_stress_time.plot(calendering_time_SN_run_1, -tau_zy_SN_run_1 * 1e-6, label=r'$\tau_{zy}$')
    #
    # # ax_force_fabric_stress_time.set_ylim(ymin=-.005, ymax=.025)
    #
    # ax_force_fabric_stress_time.legend(loc="best")
    # ax_force_fabric_stress_time.set_xlabel('Time [s]')
    # ax_force_fabric_stress_time.set_ylabel('Stress [MPa]')
    # fig_force_fabric_stress_time.tight_layout()
    #
    # fname = fig_dir + 'all_stress_time'
    # plt.savefig(fname)

    plt.show()
