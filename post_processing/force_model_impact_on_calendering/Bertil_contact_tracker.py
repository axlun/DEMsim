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

    simulation_directory = '/scratch/users/axlun/DEMsim/results/article_2/final_runs/particle_contact_model/SN_2/electrode_natural_packing_el_pl_binder_el_pl_particle/'
    particle_1, particle_2 = 3618, 4950

    time_vec, contact_data_vec = bertil_contact_tracker(simulation_directory,particle_1,particle_2)

    # =FIG 1 FORCE FOR CONTACT==========================================================================================
    fig_contact_force, ax_contact_force = plt.subplots()
    lns_total_contact = ax_contact_force.plot(time_vec,contact_data_vec[:,6], 'g', label=r'Total')
    lns_binder_contacts = ax_contact_force.plot(time_vec, contact_data_vec[:,7], 'r', label=r'Binder')
    lnd_particle_contacts = ax_contact_force.plot(time_vec,contact_data_vec[:,8],'b', label=r'Particle')
    ax_contact_force.set_ylabel("Force [N]")
    ax_contact_force.set_xlabel("time [s]")
    ax_contact_force.legend(loc='best')

    # =FIG 2 OVERLAP FOR CONTACT========================================================================================
    fig_overlap, ax_overlap = plt.subplots()
    lns_h_ = ax_overlap.plot(time_vec,contact_data_vec[:,5], 'g', label=r'h_')
    lns_h_max = ax_overlap.plot(time_vec,contact_data_vec[:,9], 'r', label=r'h_max')
    lns_b_t_ = ax_overlap.plot(time_vec,-contact_data_vec[:,20], 'b', label=r'b_t')


    ax_binder_contact = ax_overlap.twinx()
    lns_binder_contact =  ax_binder_contact.plot(time_vec,contact_data_vec[:,19], 'c+', label=r'Binder contact')
    ax_binder_contact.set_ylabel("Binder contact [-]")
    ax_overlap.set_ylabel("Overlap [m]")
    ax_overlap.set_xlabel("time [s]")
    ax_overlap.legend(loc='best')


    lns = lns_h_ + lns_h_max + lns_b_t_ + lns_binder_contact
    labs = [l.get_label() for l in lns]
    ax_overlap.legend(lns, labs, loc='best')
    # =SHOW PLOT========================================================================================================
    print('Plotting results')
    plt.show()