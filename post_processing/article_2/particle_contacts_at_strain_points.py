from Bertil_functions.Bertil_functions import bertil_data_gatherer, contact_counter_bertil
from force_model_impact_on_calendering.Bertil_mechanical_properties_multiple_runs import stress_and_linear_strain_finder

import shutil
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.style.use('axel_style')


def time_duplicate_remover(periodic_bc_data, force_fabric_tensor_data, surface_position_data):
    time, linear_strain, sxx, syy, szz, tau_xy, tau_xz, tau_yx, tau_yz, tau_zx, tau_zy \
        = stress_and_linear_strain_finder(periodic_bc_data, force_fabric_tensor_data, surface_position_data)
    rm_list = []
    for i in range(1,len(time)):
       if time[i] == time[i-1]: rm_list.append(i)
    # for i in time, linear_strain, sxx, syy, szz, tau_xy, tau_xz, tau_yx, tau_yz, tau_zx, tau_zy
    #    np.delete(i,rm_list)
    return np.delete(time,rm_list), np.delete(linear_strain,rm_list), np.delete(sxx,rm_list), np.delete(syy,rm_list), \
        np.delete(szz,rm_list), np.delete(tau_xy,rm_list), np.delete(tau_xz,rm_list), np.delete(tau_yx,rm_list), \
        np.delete(tau_yz,rm_list), np.delete(tau_zx,rm_list), np.delete(tau_zy, rm_list)


def contact_plot_processing(time, strain, contacts):
    rm_list = []
    temp_max_abs_strain = 0
    for i in range(1, len(time)):
        if temp_max_abs_strain >= abs(strain[i]):
            rm_list.append(i)
        else:
            temp_max_abs_strain = abs(strain[i])
    return np.delete(time, rm_list), np.delete(strain, rm_list), np.delete(contacts, rm_list)


class Direction:
    def __init__(self, sim_dir):
        self.time_vec, self.particle_contact_vec, self.binder_contact_vec, self.binder_particle_contact_vec \
            = contact_counter_bertil(sim_dir)

        self.force_data, self.surface_force_index, self.surface_position_index, \
            self.surface_position_data, self.periodic_BC_data,\
            self.force_fabric_tensor_data, self.kinetic_energy_data = bertil_data_gatherer(sim_dir)

        self.time, self.linear_strain, self.sxx, self.syy, self.szz, self.tau_xy, self.tau_xz, self.tau_yx,\
            self.tau_yz, self.tau_zx, self.tau_zy \
            = time_duplicate_remover(self.periodic_BC_data, self.force_fabric_tensor_data, self.surface_position_data)
        self.max_load_index = np.argmax(np.abs(self.linear_strain))

        self.time_plt, self.strain_plt, self.particle_contacts_plt = \
            contact_plot_processing(self.time, self.linear_strain, self.particle_contact_vec)

        # self.time_tension, self.linear_strain_tension, self.sxx_tension, self.syy_tension, self.szz_tension, \
        #    self.tau_xy_tension, self.tau_xz_tension, self.tau_yx_tension, self.tau_yz_tension, self.tau_zx_tension, \
        #    self.tau_zy_tension, self.time_compression, self.linear_strain_compression, self.sxx_compression, \
        #    self.syy_compression, self.szz_compression, self.tau_xy_compression, self.tau_xz, self.tau_yx, \
        #    self.tau_yz_compression, self.tau_zx_compression, self.tau_zy_compression = mech_plot_prop(sim_dir)


def contacts_at_strains_spreads_processing(sim_dir, sim_type, no_sims):
    simulation_directory = {}
    for i in range(1, no_sims + 1):
        simulation_directory[i] =\
            (Direction(sim_dir + str(i) + '/electrode_mechanical_loading_' + str(sim_type) + '_tension'),
             Direction(sim_dir + str(i) + '/electrode_mechanical_loading_' + str(sim_type) + '_compression'))
    contacts_matrix = \
        np.zeros((no_sims, np.concatenate((simulation_directory[1][0].time_plt,
                                           simulation_directory[1][1].time_plt)).size))
    for i in range(1, no_sims+1):
        contacts_matrix[i-1, :] = \
            np.concatenate((np.flip(simulation_directory[i][0].particle_contacts_plt),
                            simulation_directory[i][1].particle_contacts_plt))

    contacts_mean = np.mean(contacts_matrix, axis=0)
    contacts_std = np.std(contacts_matrix, axis=0)
    strain = \
        np.concatenate((np.flip(simulation_directory[1][0].strain_plt), simulation_directory[1][1].strain_plt))
    return strain, contacts_mean, contacts_std


class Simulation:
    def __init__(self, sim_dir, sim_type, no_sims, label='Template'):
        self.label = label
        self.stains, self.contacts_mean, self.contacts_std \
            = contacts_at_strains_spreads_processing(sim_dir, sim_type, no_sims)


if __name__ == '__main__':
    simulation_directory_SN_101 = '/scratch/users/axlun/DEMsim/results/article_2/final_runs_2/SN_101/'
    SN_101 = Simulation(simulation_directory_SN_101, 'hertz', 4, 'Reference')

    simulation_directory_SN_201 = '/scratch/users/axlun/DEMsim/results/article_2/final_runs_2/SN_201/'
    SN_201 = Simulation(simulation_directory_SN_201, 'el_pl_binder_el_pl_particle', 4, 'Material 1')

    simulation_directory_SN_301 = '/scratch/users/axlun/DEMsim/results/article_2/final_runs_2/SN_301/'
    SN_301 = Simulation(simulation_directory_SN_301, 'el_pl_binder_el_pl_particle', 4, 'Material 2')
    # ==FIGURE SAVE DIR=================================================================================================
    fig_dir = 'C:/temp/figures/particle_contacts_at_mechanical_loading/'
    try:
        shutil.rmtree(fig_dir)
        os.mkdir(fig_dir)
    except FileNotFoundError:
        print('Directory does not exist')
        os.mkdir(fig_dir)
    except:
        print('Directory could not be removed')
        quit()
    # ==================================================================================================================

    # =FIGURE CONTACTS AT STRAINS=======================================================================================
    fig_contacts, ax_contacts = plt.subplots()
    lns_SN_101 = ax_contacts.plot(SN_101.stains * 1E2, SN_101.contacts_mean, 'k-', zorder=1)
    ax_contacts.fill_between(SN_101.stains * 1E2, (SN_101.contacts_mean-SN_101.contacts_std),
                             (SN_101.contacts_mean + SN_101.contacts_std), label=SN_101.label, color='C0', alpha=1)
    lns_SN_201 = ax_contacts.plot(SN_201.stains * 1E2, SN_201.contacts_mean, 'k-', zorder=1)
    ax_contacts.fill_between(SN_201.stains * 1E2, (SN_201.contacts_mean-SN_201.contacts_std),
                             (SN_201.contacts_mean + SN_201.contacts_std), label=SN_201.label, color='C2', alpha=1)

    lns_SN_301 = ax_contacts.plot(SN_301.stains * 1E2, SN_301.contacts_mean, 'k-', zorder=1)
    ax_contacts.fill_between(SN_301.stains * 1E2, (SN_301.contacts_mean-SN_301.contacts_std),
                             (SN_301.contacts_mean + SN_301.contacts_std), label=SN_301.label, color='C3', alpha=1)
    ax_contacts.set_ylim(ymin=0)
    ax_contacts.set_ylabel('Particle contacts [-]')
    ax_contacts.set_xlabel('Strain [%]')
    ax_contacts.legend(loc='best')

    fname = fig_dir + 'contacts_at_strains'
    plt.savefig(fname)

    plt.show()
