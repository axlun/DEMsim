from Bertil_functions.Bertil_functions import *

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from os.path import exists
import shutil
import os

def local_data_gatherer(simulation_directory):
    if exists(simulation_directory + '/periodic_bc.dou'):
        periodic_BC_data = np.genfromtxt(simulation_directory + '/periodic_bc.dou', delimiter=', ')
    else:
        periodic_BC_data = np.empty([1, 1])
    force_data = np.genfromtxt(simulation_directory + '/surface_forces.dou', delimiter=', ')
    force_fabric_tensor_data = np.genfromtxt(simulation_directory + '/force_fabric_tensor.dou', delimiter=',')
    kinetic_energy_data = np.genfromtxt(simulation_directory + '/kinetic_energy.dou', delimiter=',')
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
    surface_types = [first_line[idx + 1][5:] for idx in id_idx]

    surface_indices = [i for i, surface_type in enumerate(surface_types) if surface_type == 'PointSurface']

    return force_data, surface_force_index, surface_position_index, surface_position_data, periodic_BC_data, force_fabric_tensor_data, kinetic_energy_data

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

def calendering_surface_force_func(force_data, surface_force_index):
    calendering_surface_force = force_data[:, surface_force_index[1] + 1].astype(float)
    return calendering_surface_force

def calendering_surface_pressure_func(periodic_BC_data, force_data, surface_force_index, surface_position_index):
    calendering_surface_force = force_data[:, surface_force_index[1] + 1].astype(float)
    periodic_BC_x_min = periodic_BC_data[:, 1].astype(float)
    periodic_BC_x_max = periodic_BC_data[:, 2].astype(float)
    periodic_BC_y_min = periodic_BC_data[:, 3].astype(float)
    periodic_BC_y_max = periodic_BC_data[:, 4].astype(float)
    x_side_leght = periodic_BC_x_max - periodic_BC_x_min
    y_side_leght = periodic_BC_y_max - periodic_BC_y_min
    calendering_surface_pressure = calendering_surface_force / (x_side_leght * y_side_leght)
    calendering_time = periodic_BC_data[:, 0]
    calendering_surface_position = surface_position_data[:, surface_position_index[1] + 14].astype(float)
    return calendering_time, calendering_surface_pressure, calendering_surface_position

if __name__ == '__main__':

    # ==NATUAL PACKING =================================================================================================
    # simulation_directory = '/scratch/users/axlun/DEMsim/results/final_runs/SN_run_1_q3/electrode_natural_packing_hertz'
    # simulation_directory = '/scratch/users/axlun/DEMsim/results/final_runs/SN_run_1_q_0_rmin_2_5_rmax_10/electrode_natural_packing_hertz'

    # simulation_directory = '/scratch/users/axlun/DEMsim/results/final_runs/SN_run_1/electrode_natural_packing_hertz'
    # simulation_directory = '/scratch/users/axlun/DEMsim/results/final_runs/SN_run_1_rigid_plastic_particle/electrode_natural_packing_rigid_plastic_SY_4GPa/'
    # simulation_directory = 'c:/Users/Axel/Documents/DEM/results/electrode_natural_packing_rigid_plastic_particle/SN_0/'
    simulation_directory = '/scratch/users/axlun/DEMsim/results/article_2/rigid_plastic_particle/SN_rigid_plastic_particle_N_500_dt_1e-1/electrode_natural_packing_rigid_plastic'

# ==CALENDERING=====================================================================================================
    # simulation_directory = '/scratch/users/axlun/DEMsim/results/electrode_calendering_hertz/SN_ref_run_1_5000p_btr_5_brr_15_comp_time_20_hal_105_dt_1e2_MS_1e4'
    # simulation_directory = '/scratch/users/axlun/DEMsim/results/electrode_calendering_hertz/SN_ref_run_1_10000p_btr_5_brr_15_comp_time_20_hal_105_dt_1e2_MS_1e4'
    # simulation_directory = '/scratch/users/axlun/DEMsim/results/electrode_calendering_hertz/SN_ref_run_1_5000p_btr_5_brr_15_comp_time_20_hal_105_dt_1e2_MS_1e4'
    # simulation_directory = '/scratch/users/axlun/DEMsim/results/electrode_calendering_hertz/SN_ref_run_1_5000p_btr_5_brr_15_comp_time_20_hal_105_dt_1e2_MS_1e4_rot/'
    # simulation_directory = '/scratch/users/axlun/DEMsim/results/final_runs/SN_run_1/electrode_calendering_hertz'
    # simulation_directory = '/scratch/users/axlun/DEMsim/results/final_runs/SN_run_1_q_0_rmin_3_rmax_10/electrode_calendering_hertz'
    # simulation_directory = '/scratch/users/axlun/DEMsim/results/final_runs/SN_run_1_Ebner_raw_data/electrode_calendering_hertz'
    # simulation_directory = '/scratch/users/axlun/DEMsim/results/final_runs/SN_run_1_rigid_plastic_particle/electrode_calendering_rigid_plastic/'

    # ==MECHANICAL LOADING==============================================================================================
    # simulation_directory = '/scratch/users/axlun/DEMsim/results/electrode_mechanical_loading_hertz/SN_hertz_5000p_btr_8_brr_08_dt_1e1_MS_1e2_SR_2e-3_compression'

    # ==RESTING=========================================================================================================
    # simulation_directory = '/scratch/users/axlun/DEMsim/results/electrode_resting_hertz/SN_hertz_5000p_btr_8_brr_08_dt_5e1_MS_1e4_RT_10'

    # ==PLOT PARAMETERS=================================================================================================
    fig_dir = 'c:/temp/figures/Bertil_calendering_pressure/'
    try:
        shutil.rmtree(fig_dir)
        os.mkdir(fig_dir)
    except FileNotFoundError:
        print('Directory does not exist')
        os.mkdir(fig_dir)
    except:
        print('Directory could not be removed')
        quit()
    plt.style.use('axel_style')

    # ==================================================================================================================

    force_data = []
    surface_force_index = []
    surface_position_index = []
    surface_position_data = []
    periodic_BC_data = []
    force_fabric_tensor_data = []
    kinetic_energy_data = []

    if simulation_directory.startswith("/scratch"):
        force_data, surface_force_index, surface_position_index, surface_position_data, periodic_BC_data,\
        force_fabric_tensor_data, kinetic_energy_data = bertil_data_gatherer(
            simulation_directory)
    elif simulation_directory.startswith("c:"):
        force_data, surface_force_index, surface_position_index, surface_position_data, periodic_BC_data,\
        force_fabric_tensor_data, kinetic_energy_data = local_data_gatherer(
            simulation_directory)
    else:
        print("Error with simulation directory")

    calendering_surface_force = force_data[:, surface_force_index[1] + 1].astype(float)
    bottom_surface_force = force_data[:, surface_force_index[0] + 1].astype(float)
    # periodic_BC_x_min = periodic_BC_data[:,1].astype(float)
    # periodic_BC_x_max = periodic_BC_data[:,2].astype(float)
    # periodic_BC_y_min = periodic_BC_data[:,3].astype(float)
    # periodic_BC_y_max = periodic_BC_data[:,4].astype(float)
    # x_side_length = (periodic_BC_x_max - periodic_BC_x_min)
    # y_side_length = (periodic_BC_y_max - periodic_BC_y_min)

    x_side_length = ((surface_position_data[:, surface_position_index[0] + 9]).astype(float) - (
    surface_position_data[:, surface_position_index[0] + 3]).astype(float))
    y_side_length = ((surface_position_data[:, surface_position_index[0] + 10]).astype(float) - (
    surface_position_data[:, surface_position_index[0] + 4]).astype(float))

    calendering_surface_pressure = calendering_surface_force / (x_side_length * y_side_length)
    bottom_surface_pressure = bottom_surface_force / (x_side_length * y_side_length)
    calendering_time = surface_position_data[:, -1].astype(float)
    calendering_surface_position = surface_position_data[:, surface_position_index[1] + 14].astype(float)

    # =================================STRESS IN Z-DIRECTION============================================================
    vol = x_side_length * y_side_length * calendering_surface_position
    sxx = -force_fabric_tensor_data[:, 1] / vol
    syy = -force_fabric_tensor_data[:, 5] / vol
    szz = -force_fabric_tensor_data[:, 9] / vol
    tau_xy = -force_fabric_tensor_data[:, 2] / vol
    tau_xz = -force_fabric_tensor_data[:, 3] / vol

    tau_yx = -force_fabric_tensor_data[:, 4] / vol
    tau_yz = -force_fabric_tensor_data[:, 6] / vol

    tau_zx = -force_fabric_tensor_data[:, 7] / vol
    tau_zy = -force_fabric_tensor_data[:, 8] / vol


    # ===FIG 1 CALENDERING SURFACE PRESSURE SURFACE POSITION TIME=======================================================
    fig_calendering_surface_pressure, ax_calendering_surface_pressure = plt.subplots()
    ax_calendering_surface_pressure.set_ylabel("Calendering surface pressure [MPa]")
    ax_calendering_surface_pressure.set_xlabel("Time [s]")
    lns_calendering_surface_pressure = ax_calendering_surface_pressure.plot(calendering_time,
                                                                            calendering_surface_pressure * 1E-6,
                                                                            label='Pressure')
    ax_calendering_surface_position = ax_calendering_surface_pressure.twinx()
    ax_calendering_surface_position.set_ylabel("Calendering surface position [m]")
    lns_calendering_surface_position = ax_calendering_surface_position.plot(calendering_time,
                                                                            calendering_surface_position, 'c',
                                                                            label='Position')
    # ax_calendering_surface_pressure.set_title('Calendering surface pressure')

    lns = lns_calendering_surface_pressure + lns_calendering_surface_position
    labs = [l.get_label() for l in lns]
    ax_calendering_surface_pressure.legend(lns, labs, loc='best')
    fig_calendering_surface_pressure.tight_layout()

    fname = fig_dir + 'calendering_surface_pressure_surface_position_time'
    plt.savefig(fname)

    # ===FIG 2 CALENDERING SURFACE PRESSURE POSITION====================================================================
    calendering_initiate_index, calendering_break_index = calendering_break_index_func(calendering_surface_position)
    fig_calendering_process, ax_calendering_process = plt.subplots()
    ax_calendering_process.plot(calendering_surface_position[calendering_initiate_index:(calendering_break_index - 1)],
                                calendering_surface_pressure[
                                calendering_initiate_index:(calendering_break_index - 1)] * 1e-6)
    ax_calendering_process.set_ylabel("Calendering surface pressure [MPa]")
    ax_calendering_process.set_xlabel("Calendering surface position [m]")

    # ax_calendering_process.set_title('Calendering surface pressure')
    fig_calendering_process.tight_layout()

    fname = fig_dir + 'calendering_surface_pressure_surface_position'
    plt.savefig(fname)

    # ===FIG 3 CALENDERING SURFACE PRESSURE SURFACE POSITION FULL=======================================================
    fig_full_calendering_sim, ax_full_calendering_sim = plt.subplots()
    ax_full_calendering_sim.plot(calendering_surface_position[:], calendering_surface_pressure[:] * 1E-6)
    ax_full_calendering_sim.set_ylabel("Calendering surface pressure [MPa]")
    ax_full_calendering_sim.set_xlabel("Calendering surface position [m]")
    #    ax_full_calendering_sim.set_title('calendering surface pressure')
    fig_full_calendering_sim.tight_layout()

    fname = fig_dir + 'calendering_surface_pressure_surface_position_full'
    plt.savefig(fname)

    # ===FIG 4 STRESS ZZ SURFACE POSITION===============================================================================

    fig_sigma_zz_calendering_surface_position, ax_sigma_zz_calendering_surface_position = plt.subplots()
    ax_sigma_zz_calendering_surface_position.plot(calendering_surface_position[:], szz)
    ax_sigma_zz_calendering_surface_position.set_ylabel("Stress in z [Pa]")
    ax_sigma_zz_calendering_surface_position.set_xlabel("Calendering surface position [m]")
    # ax_sigma_zz_calendering_surface_position.set_title("Stress in Z to calendering surface position")
    fig_sigma_zz_calendering_surface_position.tight_layout()

    # ===FIG 5 STRESS ZZ SURFACE PRESSURE SURFACE POSITION TIME=========================================================
    fig_sigma_zz_time, ax_sigma_zz_time = plt.subplots()
    lns_sigma_zz_time = ax_sigma_zz_time.plot(calendering_time, -szz * 1e-6, label=r'$\sigma_{zz}$')
    lsn_cal_surf_pre_time = ax_sigma_zz_time.plot(calendering_time, calendering_surface_pressure * 1e-6,
                                                  label=r'Surface pressure')
    ax_sigma_zz_time.set_ylabel("Stress in z [MPa]")
    ax_sigma_zz_time.set_xlabel("Time [s]")
    ax_calendering_surface_position_time = ax_sigma_zz_time.twinx()
    ax_calendering_surface_position_time.set_ylabel("Calendering surface position [m]")
    lns_calendering_surface_position_time = ax_calendering_surface_position_time.plot(calendering_time,
                                                                                      calendering_surface_position, 'c',
                                                                                      label='Surface position')
    lns2 = lns_sigma_zz_time + lsn_cal_surf_pre_time + lns_calendering_surface_position_time
    labs2 = [l.get_label() for l in lns2]
    ax_sigma_zz_time.legend(lns2, labs2, loc="best")
    fig_sigma_zz_time.tight_layout()

    # ===FIG 6 STRESS ZZ SURFACE PRESSURE TIME==========================================================================
    fig_calendering_surface_pressure_sigma_zz, ax_sigma_zz_time2 = plt.subplots()
    lns5 = ax_sigma_zz_time2.plot(calendering_time, -szz, label=r'$\sigma_{zz}$')
    ax_sigma_zz_time2.set_ylabel("Stress [Pa]")
    ax_sigma_zz_time2.set_xlabel("Time [s]")
    lns6 = ax_sigma_zz_time2.plot(calendering_time, calendering_surface_pressure[:], label='Surface pressure')
    #  ax9.set_title('calendering surface pressure')
    lns3 = lns5 + lns6
    labs3 = [l.get_label() for l in lns3]
    ax_sigma_zz_time2.legend(lns3, labs3, loc="best")
    fig_calendering_surface_pressure_sigma_zz.tight_layout()

    # ===FIG 7 ALL SURFACE PRESSURES STRESS ZZ SURFACE POSITION TIME====================================================
    fig_calendering_surface_bottom_surface, ax_calendering_surface_bottom_surface = plt.subplots()
    lns_calendering_surface_pressure_2 = ax_calendering_surface_bottom_surface.plot(calendering_time,
                                                                                    calendering_surface_pressure * 1e-6,
                                                                                    'g', label=r'Surface pressure')
    lns_bottom_surface_pressure = ax_calendering_surface_bottom_surface.plot(calendering_time,
                                                                             bottom_surface_pressure * 1e-6,
                                                                             label=r'Bottom pressure')
    lns_macroscopic_stress = ax_calendering_surface_bottom_surface.plot(calendering_time, -szz * 1e-6,
                                                                        label=r'$\sigma_{zz}$')

    ax_calendering_surface_position_2 = ax_calendering_surface_bottom_surface.twinx()
    ax_calendering_surface_position_2.set_ylabel("Calendering surface position [m]")
    #    ax_calendering_surface_position_2.set_ylim([0,200])

    lns_calendering_surface_position_2 = ax_calendering_surface_position_2.plot(calendering_time,
                                                                                calendering_surface_position, 'c',
                                                                                label='Surface position')
    ax_calendering_surface_bottom_surface.set_ylabel("Stress in z [MPa]")
    ax_calendering_surface_bottom_surface.set_xlabel("Time [s]")
    # ax_calendering_surface_bottom_surface.set_title("Pressure on calendering surface and bottom surface")
    lns_sum = lns_calendering_surface_pressure_2 + lns_bottom_surface_pressure + lns_macroscopic_stress + lns_calendering_surface_position_2
    labs4 = [l.get_label() for l in lns_sum]
    # ax_calendering_surface_bottom_surface.set_ylim(ymin=-.02,ymax=.12)
    ax_calendering_surface_bottom_surface.legend(lns_sum, labs4, loc="best")
    fig_calendering_surface_bottom_surface.tight_layout()

    fname = fig_dir + 'all_surface_pressure_stess_ZZ_surface_position_time'
    plt.savefig(fname)

    # ===FIG 8 KINETIC  ENERGY===============================================================================================
    fig_KE, ax_KE = plt.subplots()
    lns5 = ax_KE.plot(kinetic_energy_data[:, -1], kinetic_energy_data[:, 0], label=r'KE')
    ax_KE.set_ylabel("Kinetic energy [J]")
    ax_KE.set_xlabel("Time [s]")
    # ax_KE.set_title("Kinetic energy of particles")
    ax_KE.legend(loc="best")
    fig_KE.tight_layout()

    fname = fig_dir + 'KE'
    plt.savefig(fname)

    # ===FIG 9 STRESS XX YY ZZ TIME===============================================================================================
    fig_force_fabric_stress_time, ax_force_fabric_stress_time = plt.subplots()
    lns_sxx = ax_force_fabric_stress_time.plot(calendering_time, -sxx * 1e-6, label=r'$\sigma_{xx}$')
    lns_syy = ax_force_fabric_stress_time.plot(calendering_time, -syy * 1e-6, label=r'$\sigma_{yy}$')
    lns_szz = ax_force_fabric_stress_time.plot(calendering_time, -szz * 1e-6, label=r'$\sigma_{zz}$')
    lns_tau_xy = ax_force_fabric_stress_time.plot(calendering_time, -tau_xy * 1e-6, label=r'$\tau_{xy}$')
    lns_tau_xz = ax_force_fabric_stress_time.plot(calendering_time, -tau_xz * 1e-6, label=r'$\tau_{xz}$')
    lns_tau_yx = ax_force_fabric_stress_time.plot(calendering_time, -tau_yx * 1e-6, label=r'$\tau_{yx}$')
    lns_tau_yz = ax_force_fabric_stress_time.plot(calendering_time, -tau_yz * 1e-6, label=r'$\tau_{yz}$')
    lns_tau_zx = ax_force_fabric_stress_time.plot(calendering_time, -tau_zx * 1e-6, label=r'$\tau_{zx}$')
    lns_tau_zy = ax_force_fabric_stress_time.plot(calendering_time, -tau_zy * 1e-6, label=r'$\tau_{zy}$')

    # ax_force_fabric_stress_time.set_ylim(ymin=-.005, ymax=.025)

    ax_force_fabric_stress_time.legend(loc="best")
    ax_force_fabric_stress_time.set_xlabel('Time [s]')
    ax_force_fabric_stress_time.set_ylabel('Stress [MPa]')
    fig_force_fabric_stress_time.tight_layout()

    fname = fig_dir + 'all_stress_time'
    plt.savefig(fname)

    plt.show()
