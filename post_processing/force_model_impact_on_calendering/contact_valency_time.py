from Bertil_functions.Bertil_functions import *
from calendering_pressure import local_data_gatherer, bertil_data_gatherer

import numpy as np

import matplotlib.pyplot as plt
import re
import matplotlib

matplotlib.style.use('classic')



if __name__ == '__main__':

#    simulation_directory = '/scratch/users/axlun/DEMsim/results/electrode_calendaring_restart/SN_hertz_5000p_btr_00_hal_6'
    simulation_directory = '/scratch/users/axlun/DEMsim/results/electrode_calendaring/SN_hertz_3000p_btr_065'#_new_Ft_b'
#    simulation_directory = '/scratch/users/axlun/DEMsim/results/electrode_calendaring/SN00_400p_plastic_binder'

    time_vec,particle_contact_valency_mat,binder_contact_valency_mat = contact_valency_bertil(simulation_directory)

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

    for i in range(0,len(particle_contact_valency_mat[:])):
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

    fig_particle_contact_valency, ax_particle_contact_valency = plt.subplots()
    ax_particle_contact_valency.set_ylabel('Number of particles [-]')
    ax_particle_contact_valency.set_xlabel('Time [s]')

    lns_particle_contact_valency_0 = ax_particle_contact_valency.plot(time_vec[:], particle_contact_valency_0_vec, lw=1,
                                             label=r'0')
    lns_particle_contact_valency_1 = ax_particle_contact_valency.plot(time_vec[:], particle_contact_valency_1_vec, lw=1,
                                                label=r'1')
    lns_particle_contact_valency_2 = ax_particle_contact_valency.plot(time_vec[:], particle_contact_valency_2_vec, lw=1,
                                                label=r'2')
    lns_particle_contact_valency_3 = ax_particle_contact_valency.plot(time_vec[:], particle_contact_valency_3_vec, lw=1,
                                                label=r'3')
    lns_particle_contact_valency_4 = ax_particle_contact_valency.plot(time_vec[:], particle_contact_valency_4_vec, lw=1,
                                                label=r'4')
    lns_particle_contact_valency_5 = ax_particle_contact_valency.plot(time_vec[:], particle_contact_valency_5_vec, lw=1,
                                                label=r'5')
    lns_particle_contact_valency_6 = ax_particle_contact_valency.plot(time_vec[:], particle_contact_valency_6_vec, lw=1,
                                                label=r'6')
    lns_particle_contact_valency_7 = ax_particle_contact_valency.plot(time_vec[:], particle_contact_valency_7_vec, lw=1,
                                                label=r'7')

    ax_particle_contact_valency.set_title('Contact valency of particle particle contacts')
    plt.legend(loc='best')

    fig_binder_contact_valency, ax_binder_contact_valency = plt.subplots()
    ax_binder_contact_valency.set_ylabel('Number of paticles [-]')
    ax_binder_contact_valency.set_xlabel('Time [s]')

    lns_binder_contact_valency_0 = ax_binder_contact_valency.plot(time_vec[:], binder_contact_valency_0_vec, lw=1,
                                             label=r'0')
    lns_binder_contact_valency_1 = ax_binder_contact_valency.plot(time_vec[:], binder_contact_valency_1_vec, lw=1,
                                                label=r'1')
    lns_binder_contact_valency_2 = ax_binder_contact_valency.plot(time_vec[:], binder_contact_valency_2_vec, lw=1,
                                                label=r'2')
    lns_binder_contact_valency_3 = ax_binder_contact_valency.plot(time_vec[:], binder_contact_valency_3_vec, lw=1,
                                                label=r'3')
    lns_binder_contact_valency_4 = ax_binder_contact_valency.plot(time_vec[:], binder_contact_valency_4_vec, lw=1,
                                                label=r'4')
    lns_binder_contact_valency_5 = ax_binder_contact_valency.plot(time_vec[:], binder_contact_valency_5_vec, lw=1,
                                                label=r'5')
    lns_binder_contact_valency_6 = ax_binder_contact_valency.plot(time_vec[:], binder_contact_valency_6_vec, lw=1,
                                                              label=r'6')
    lns_binder_contact_valency_7 = ax_binder_contact_valency.plot(time_vec[:], binder_contact_valency_7_vec, lw=1,
                                                              label=r'7')

    ax_binder_contact_valency.set_title('Contact valency of binder contacts')
    plt.legend(loc='best')

    plt.show()


