from Bertil_functions.Bertil_functions import *
from Bertil_calendering_pressure import local_data_gatherer, bertil_data_gatherer

import numpy as np

import matplotlib.pyplot as plt
import re
import matplotlib
import os
import shutil

if __name__ == '__main__':

    # simulation_directory = '/scratch/users/axlun/DEMsim/results/article_4/particle_fracture/2/electrode_calendering_binder_fracture_005/'
    # simulation_directory = '/scratch/users/axlun/DEMsim/results/article_4/particle_fracture/2/electrode_calendering_binder_fracture_01/'
    # simulation_directory = '/scratch/users/axlun/DEMsim/results/article_4/particle_fracture/2/electrode_calendering_binder_fracture_02/'

    # simulation_directory = '/scratch/users/axlun/DEMsim/results/article_4/particle_fracture/2/electrode_calendering_only_binder_fracture_005/'
    # simulation_directory = '/scratch/users/axlun/DEMsim/results/article_4/particle_fracture/2/electrode_calendering_only_binder_fracture_01/'
    # simulation_directory = '/scratch/users/axlun/DEMsim/results/article_4/particle_fracture/2/electrode_calendering_only_binder_fracture_02/'

    # simulation_directory = '/scratch/users/axlun/DEMsim/results/article_4/particle_fracture/2/electrode_mechanical_loading_only_binder_fracture_005'
    simulation_directory = '/scratch/users/axlun/DEMsim/results/article_4/particle_fracture/2/electrode_mechanical_loading_only_binder_fracture_01_tension'
    # simulation_directory = '/scratch/users/axlun/DEMsim/results/article_4/particle_fracture/2/electrode_mechanical_loading_only_binder_fracture_02_compression'
    # simulation_directory = '/scratch/users/axlun/DEMsim/results/article_4/particle_fracture/2/electrode_mechanical_loading_only_binder_fracture_02'

    # simulation_directory = '/scratch/users/axlun/DEMsim/results/article_4/particle_fracture/2/electrode_mechanical_loading_tension_cpy'

    fig_dir = 'C:/temp/figures/Bertil_fractured_binder_contacts/'
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

    time_vec, binder_contact_vec,active_binder_contact_vec, fractured_binder_contact_vec = fractured_binder_counter(simulation_directory)

#===FIG 1 NUMBER OF CONTACTS=============================================================================================
    fig_contact_distribution, ax_contact_distribution = plt.subplots()
    lns_binder_contacts = ax_contact_distribution.plot(time_vec, binder_contact_vec, label=r'Binder contacts')
    lns_active_binder_contacts = ax_contact_distribution.plot(time_vec, active_binder_contact_vec, label=r'Active binder contacts')
    lns_fractured_binder_contacts = ax_contact_distribution.plot(time_vec, fractured_binder_contact_vec, label=r'Fractured binder contacts')
    lns_active_fractured_binder_contacts = ax_contact_distribution.plot(time_vec, active_binder_contact_vec + fractured_binder_contact_vec, label=r'Active + fractured binder contacts')
    ax_contact_distribution.set_ylabel("Number of contacts [-]")
    ax_contact_distribution.set_xlabel("Time [s]")
    ax_contact_distribution.legend(loc='best')
    ax_contact_distribution.set_ylim(ymin=0)
    ax_contact_distribution.set_xlim(xmin=15, xmax=35)
    fig_contact_distribution.tight_layout()

    fig_normalised_fractured_contacts, ax_normalised_fractured_contacts = plt.subplots()
    lns_normalised_fractured_contacts = ax_normalised_fractured_contacts.plot(time_vec, fractured_binder_contact_vec/binder_contact_vec)
    ax_normalised_fractured_contacts.set_ylabel("Fractured binder contacts over total [-]")
    ax_normalised_fractured_contacts.set_xlabel("Time [s]")
    ax_normalised_fractured_contacts.set_xlim(xmin=15, xmax=35)
    ax_normalised_fractured_contacts.set_ylim(ymin=0)
    fig_normalised_fractured_contacts.tight_layout()


    fname = fig_dir + 'number_of_contacts'
    plt.savefig(fname)

    print('Plotting results')
    plt.show()