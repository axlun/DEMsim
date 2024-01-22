from Bertil_functions.Bertil_functions import *
from Bertil_calendering_pressure import local_data_gatherer, bertil_data_gatherer

import numpy as np

import matplotlib.pyplot as plt
import sys
import os
import re
import matplotlib
import pandas as pd


if __name__ == '__main__':
    plt.style.use('axel_style')

    # simulation_directory = '/scratch/users/axlun/DEMsim/results/final_runs/SN_run_1/electrode_mechanical_loading_hertz_compression'
    # simulation_directory = '/scratch/users/axlun/DEMsim/results/final_runs/SN_run_1/electrode_calendering_hertz'
    # simulation_directory = '/scratch/users/axlun/DEMsim/results/article_2/viscoelastic_testing/SN_1/electrode_relaxation_el_pl_binder_el_pl_particle_compression_01/'
    # simulation_directory = '/scratch/users/axlun/DEMsim/results/article_2/viscoelastic_testing/SN_1/electrode_calendering_el_pl_binder_el_pl_particle/'
    # simulation_directory = '/scratch/users/axlun/DEMsim/results/article_2/viscoelastic_testing/SN_4/electrode_natural_packing_el_pl_binder_el_pl_particle/'
    # simulation_directory = '/scratch/users/axlun/DEMsim/results/article_2/viscoelastic_testing/SN_1/electrode_natural_packing_el_pl_binder_el_pl_particle/'
    # particle_1 = 857

    # simulation_directory = '/scratch/users/axlun/DEMsim/results/article_2/final_runs/particle_contact_model/SN_2/electrode_natural_packing_el_pl_binder_el_pl_particle/'
    # particle_1 = 4950

    # simulation_directory = '/scratch/users/axlun/DEMsim/results/article_2/final_runs/particle_contact_model/SN_1/electrode_calendering_el_pl_binder_el_pl_particle/'
    # simulation_directory = '/scratch/users/axlun/DEMsim/results/article_2/final_runs/particle_contact_model/SN_1/electrode_natural_packing_el_pl_binder_el_pl_particle/'
    # particle_1 = 3676
    # particle_1 = 4003


    # simulation_directory = '/scratch/users/axlun/DEMsim/results/article_2/final_runs/particle_contact_model/SN_1/electrode_calendering_el_pl_binder_el_pl_particle/'
    # particle_1 = 4273

    simulation_directory = '/scratch/users/axlun/DEMsim/results/article_2/final_runs/particle_contact_model/SN_1_2_2/electrode_mechanical_loading_el_pl_binder_el_pl_particle_compression/'
    particle_1 = 1581

    time_vec, contact_data_vec = bertil_single_contact_tracker(simulation_directory,particle_1)
    print(time_vec)
    print(contact_data_vec)
    # # =FIG 1 FORCE FOR CONTACT==========================================================================================
    fig_contacts_time, ax_contacts_time = plt.subplots()
    lns_total_contact = ax_contacts_time.plot(time_vec,contact_data_vec[:,0], 'g', label=r'Total')
    lns_binder_contacts = ax_contacts_time.plot(time_vec, contact_data_vec[:,1], 'r', label=r'Binder')
    # lnd_particle_contacts = ax_contacts_time.plot(time_vec,contact_data_vec[:,8],'b', label=r'Particle')
    ax_contacts_time.set_ylabel("Number of contacts [-]")
    ax_contacts_time.set_xlabel("time [s]")
    ax_contacts_time.legend(loc='best')
    #
    # # =FIG 2 OVERLAP FOR CONTACT========================================================================================
    # fig_overlap, ax_overlap = plt.subplots()
    # lns_h_ = ax_overlap.plot(time_vec,contact_data_vec[:,5], 'g', label=r'h_')
    # lns_h_max = ax_overlap.plot(time_vec,contact_data_vec[:,9], 'r', label=r'h_max')
    # lns_b_t_ = ax_overlap.plot(time_vec,-contact_data_vec[:,20], 'b', label=r'b_t')
    #
    #
    # ax_binder_contact = ax_overlap.twinx()
    # lns_binder_contact =  ax_binder_contact.plot(time_vec,contact_data_vec[:,19], 'c+', label=r'Binder contact')
    # ax_binder_contact.set_ylabel("Binder contact [-]")
    # ax_overlap.set_ylabel("Overlap [m]")
    # ax_overlap.set_xlabel("time [s]")
    # ax_overlap.legend(loc='best')


    # lns = lns_h_ + lns_h_max + lns_b_t_ + lns_binder_contact
    # labs = [l.get_label() for l in lns]
    # ax_overlap.legend(lns, labs, loc='best')
    # =SHOW PLOT========================================================================================================
    plt.show()