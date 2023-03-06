from Bertil_functions.Bertil_functions import *
from Bertil_calendering_pressure import local_data_gatherer, bertil_data_gatherer

import numpy as np

import matplotlib.pyplot as plt
import re
import matplotlib
import os
import shutil

if __name__ == '__main__':

    # ==NATURAL PACKING=================================================================================================
    # simulation_directory = '/scratch/users/axlun/DEMsim/results/electrode_natural_packing_hertz/SN_hertz_5000p_btr_5_brr_05_dt_1e0_MS_1e0'

    # ==CALENDERING=====================================================================================================
    # simulation_directory = '/scratch/users/axlun/DEMsim/results/electrode_calendering_hertz/SN_ref_run_1_5000p_btr_5_brr_8_comp_time_20_hal_105_dt_1e2_MS_1e4'
    # simulation_directory = '/scratch/users/axlun/DEMsim/results/final_runs/SN_run_1_El_Pl/electrode_calendering_hertz'
    simulation_directory = '/scratch/users/axlun/DEMsim/results/final_runs/SN_run_1_rigid_plastic_particle/electrode_calendering_rigid_plastic/'


    # ==MECHANICAL LOADING==============================================================================================
    # simulation_directory = '/scratch/users/axlun/DEMsim/results/electrode_mechanical_loading_hertz/SN_hertz_5000p_btr_5_brr_05_dt_5e1_MS_1e2_SR_2e-3_compression'
    # simulation_directory = '/scratch/users/axlun/DEMsim/results/final_runs/SN_run_1_El_Pl/electrode_mechanical_loading_hertz_compression'

    time_vec, particle_contact_valency_mat, binder_contact_valency_mat = contact_valency_bertil(simulation_directory)

    # ==PLOT PARAMETERS=================================================================================================
    fig_dir = 'C:/temp/figures/Bertil_contact_valency/'
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

    particle_contact_valency_0_vec = []
    particle_contact_valency_1_vec = []
    particle_contact_valency_2_vec = []
    particle_contact_valency_3_vec = []
    particle_contact_valency_4_vec = []
    particle_contact_valency_5_vec = []
    particle_contact_valency_6_vec = []
    particle_contact_valency_7_vec = []

    binder_contact_valency_0_vec = []
    binder_contact_valency_1_vec = []
    binder_contact_valency_2_vec = []
    binder_contact_valency_3_vec = []
    binder_contact_valency_4_vec = []
    binder_contact_valency_5_vec = []
    binder_contact_valency_6_vec = []
    binder_contact_valency_7_vec = []

    for i in range(0, len(particle_contact_valency_mat[:])):
        particle_contact_valency_0_vec.append(int(particle_contact_valency_mat[i][0]))
        particle_contact_valency_1_vec.append(int(particle_contact_valency_mat[i][1]))
        particle_contact_valency_2_vec.append(int(particle_contact_valency_mat[i][2]))
        particle_contact_valency_3_vec.append(int(particle_contact_valency_mat[i][3]))
        particle_contact_valency_4_vec.append(int(particle_contact_valency_mat[i][4]))
        particle_contact_valency_5_vec.append(int(particle_contact_valency_mat[i][5]))
        particle_contact_valency_6_vec.append(int(particle_contact_valency_mat[i][6]))
        particle_contact_valency_7_vec.append(int(particle_contact_valency_mat[i][7]))

        binder_contact_valency_0_vec.append(int(binder_contact_valency_mat[i][0]))
        binder_contact_valency_1_vec.append(int(binder_contact_valency_mat[i][1]))
        binder_contact_valency_2_vec.append(int(binder_contact_valency_mat[i][2]))
        binder_contact_valency_3_vec.append(int(binder_contact_valency_mat[i][3]))
        binder_contact_valency_4_vec.append(int(binder_contact_valency_mat[i][4]))
        binder_contact_valency_5_vec.append(int(binder_contact_valency_mat[i][5]))
        binder_contact_valency_6_vec.append(int(binder_contact_valency_mat[i][6]))
        binder_contact_valency_7_vec.append(int(binder_contact_valency_mat[i][7]))

    # ==PARTICLE CONTACT VALENCY========================================================================================
    fig_particle_contact_valency, ax_particle_contact_valency = plt.subplots()
    ax_particle_contact_valency.set_ylabel('Number of particles [-]')
    ax_particle_contact_valency.set_xlabel('Time [s]')

    lns_particle_contact_valency_0 = ax_particle_contact_valency.plot(time_vec[:], particle_contact_valency_0_vec,
                                                                      label=r'0')
    lns_particle_contact_valency_1 = ax_particle_contact_valency.plot(time_vec[:], particle_contact_valency_1_vec,
                                                                      label=r'1')
    lns_particle_contact_valency_2 = ax_particle_contact_valency.plot(time_vec[:], particle_contact_valency_2_vec,
                                                                      label=r'2')
    lns_particle_contact_valency_3 = ax_particle_contact_valency.plot(time_vec[:], particle_contact_valency_3_vec,
                                                                      label=r'3')
    lns_particle_contact_valency_4 = ax_particle_contact_valency.plot(time_vec[:], particle_contact_valency_4_vec,
                                                                      label=r'4')
    lns_particle_contact_valency_5 = ax_particle_contact_valency.plot(time_vec[:], particle_contact_valency_5_vec,
                                                                      label=r'5')
    lns_particle_contact_valency_6 = ax_particle_contact_valency.plot(time_vec[:], particle_contact_valency_6_vec,
                                                                      label=r'6')
    lns_particle_contact_valency_7 = ax_particle_contact_valency.plot(time_vec[:], particle_contact_valency_7_vec,
                                                                      label=r'>7')

    # ax_particle_contact_valency.set_title('Contact valency of particle particle contacts')
    ax_particle_contact_valency.set_ylim(ymin=0)
    plt.legend(loc='upper center', ncol=8, bbox_to_anchor=(.5, 1.2), title='Number of particle contacts')
    fig_particle_contact_valency.tight_layout()

    fname = fig_dir + 'particle_contact_valency'
    plt.savefig(fname)

    # ==BINDER CONTACT VALENCY==========================================================================================
    fig_binder_contact_valency, ax_binder_contact_valency = plt.subplots()
    ax_binder_contact_valency.set_ylabel('Number of particles [-]')
    ax_binder_contact_valency.set_xlabel('Time [s]')

    lns_binder_contact_valency_0 = ax_binder_contact_valency.plot(time_vec[:], binder_contact_valency_0_vec,
                                                                  label=r'0')
    lns_binder_contact_valency_1 = ax_binder_contact_valency.plot(time_vec[:], binder_contact_valency_1_vec,
                                                                  label=r'1')
    lns_binder_contact_valency_2 = ax_binder_contact_valency.plot(time_vec[:], binder_contact_valency_2_vec,
                                                                  label=r'2')
    lns_binder_contact_valency_3 = ax_binder_contact_valency.plot(time_vec[:], binder_contact_valency_3_vec,
                                                                  label=r'3')
    lns_binder_contact_valency_4 = ax_binder_contact_valency.plot(time_vec[:], binder_contact_valency_4_vec,
                                                                  label=r'4')
    lns_binder_contact_valency_5 = ax_binder_contact_valency.plot(time_vec[:], binder_contact_valency_5_vec,
                                                                  label=r'5')
    lns_binder_contact_valency_6 = ax_binder_contact_valency.plot(time_vec[:], binder_contact_valency_6_vec,
                                                                  label=r'6')
    lns_binder_contact_valency_7 = ax_binder_contact_valency.plot(time_vec[:], binder_contact_valency_7_vec,
                                                                  label=r'>7')

    # ax_binder_contact_valency.set_title('Contact valency of binder contacts')
    ax_binder_contact_valency.set_ylim(ymin=0)

    plt.legend(loc='upper center', ncol=8, bbox_to_anchor=(.5, 1.2), title='Number of binder contacts')
    fig_binder_contact_valency.tight_layout()

    fname = fig_dir + 'binder_contact_valency'
    plt.savefig(fname)

    # ==SHOWING PLOT====================================================================================================
    print('Plotting')
    plt.show()
