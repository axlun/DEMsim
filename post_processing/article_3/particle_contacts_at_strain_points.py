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
            (Direction(sim_dir + str(i) + '/swelling_electrode_mechanical_loading_' + str(sim_type) + 'tension'),
             Direction(sim_dir + str(i) + '/swelling_electrode_mechanical_loading_' + str(sim_type) + 'compression'))
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
        self.strains, self.contacts_mean, self.contacts_std \
            = contacts_at_strains_spreads_processing(sim_dir, sim_type, no_sims)


if __name__ == '__main__':

    no_sims = 4 #TODO:: run with 4 sims

    SOC_50 = Simulation('/scratch/users/axlun/DEMsim/results/article_3/final_runs/',
                        'ss_1.03228_material_scaling_', no_sims, '50')

    SOC_0 = Simulation('/scratch/users/axlun/DEMsim/results/article_3/final_runs/', '', no_sims, '0')

    SOC_100 = Simulation('/scratch/users/axlun/DEMsim/results/article_3/final_runs/',
                         'ss_1.06266_material_scaling_', no_sims, '100')


    # swelling = Simulation('/scratch/users/axlun/DEMsim/results/article_3/final_runs/',
    #                       'ss_1.06266_', no_sims , 'swelling')

    # material_deg = Simulation('/scratch/users/axlun/DEMsim/results/article_3/final_runs/',
    #                           'material_scaling_', no_sims , 'material degradation')

    cycle_1 = Simulation('/scratch/users/axlun/DEMsim/results/article_3/final_runs/',
                         'cycle_1_', no_sims, '1')

    cycle_3 = Simulation('/scratch/users/axlun/DEMsim/results/article_3/final_runs/',
                         'cycle_3_', no_sims, '3')

    cycle_10 = Simulation('/scratch/users/axlun/DEMsim/results/article_3/final_runs/',
                          'cycle_10_', no_sims, '10')

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

    # =FIGURE SOC=======================================================================================================
    fig_soc, ax_soc = plt.subplots()

    lns_soc_100 = ax_soc.plot(SOC_100.strains * 1E2, SOC_100.contacts_mean, 'k-', zorder=3)
    ax_soc.fill_between(SOC_100.strains * 1E2,
                        (SOC_100.contacts_mean-SOC_100.contacts_std),
                        (SOC_100.contacts_mean + SOC_100.contacts_std),
                        label=SOC_100.label, color='C0', alpha=1)

    lns_soc_50 = ax_soc.plot(SOC_50.strains * 1E2, SOC_50.contacts_mean, 'k-', zorder=2)
    ax_soc.fill_between(SOC_50.strains * 1E2,
                        (SOC_50.contacts_mean-SOC_50.contacts_std),
                        (SOC_50.contacts_mean + SOC_50.contacts_std),
                        label=SOC_50.label, color='C2', alpha=1)

    lns_soc_0 = ax_soc.plot(SOC_0.strains * 1E2, SOC_0.contacts_mean, 'k-', zorder=1)
    ax_soc.fill_between(SOC_0.strains * 1E2,
                        (SOC_0.contacts_mean-SOC_0.contacts_std),
                        (SOC_0.contacts_mean + SOC_0.contacts_std),
                        label=SOC_0.label, color='C1', alpha=1)

    ax_soc.set_ylim(ymin=0)
    ax_soc.set_xlim(xmin=-2.2, xmax=2.2)
    ax_soc.set_ylabel('Particle contacts [-]')
    ax_soc.set_xlabel('Strain [%]')
    ax_soc.legend(loc='best', title='SOC [%]')
    fig_soc.tight_layout()
    fname = fig_dir + 'SOC'
    plt.savefig(fname)

    """
    # =FIGURE SPLIT=====================================================================================================
    fig_split, ax_split= plt.subplots()


    ax_split.plot(SOC_100.strains * 1E2, SOC_100.contacts_mean, 'k-', zorder=1)
    ax_split.fill_between(SOC_100.strains * 1E2, (SOC_100.contacts_mean-SOC_100.contacts_std),
                          (SOC_100.contacts_mean + SOC_100.contacts_std),
                          label='Swelling + material degradation', color='C2', alpha=1)

    ax_split.plot(swelling.strains * 1E2, swelling.contacts_mean, 'k-', zorder=2)
    ax_split.fill_between(swelling.strains * 1E2, (swelling.contacts_mean-swelling.contacts_std),
                          (swelling.contacts_mean + swelling.contacts_std), label='Swelling', color='C1', alpha=1)

    ax_split.plot(SOC_0.strains * 1E2, SOC_0.contacts_mean, 'k-', zorder=1)
    ax_split.fill_between(SOC_0.strains * 1E2, (SOC_0.contacts_mean-SOC_0.contacts_std),
                          (SOC_0.contacts_mean + SOC_0.contacts_std), label='Reference', color='C0', alpha=1)


    ax_split.set_ylim(ymin=0)
    ax_split.set_ylabel('Particle contacts [-]')
    ax_split.set_xlabel('Strain [%]')
    ax_split.legend(loc='best')

    fname = fig_dir + 'split'
    plt.savefig(fname)
    """

    # =FIGURE CYCLING===================================================================================================
    fig_cycling, ax_cycling= plt.subplots()

    lns_cycle_10 = ax_cycling.plot(cycle_10.strains * 1E2, cycle_10.contacts_mean, 'k-', zorder=3)
    ax_cycling.fill_between(cycle_10.strains * 1E2,
                            (cycle_10.contacts_mean-cycle_10.contacts_std),
                            (cycle_10.contacts_mean + cycle_10.contacts_std),
                            label=cycle_10.label, color='C6', alpha=1)

    lns_cycle_3 = ax_cycling.plot(cycle_3.strains * 1E2, cycle_3.contacts_mean, 'k-', zorder=3)
    ax_cycling.fill_between(cycle_3.strains * 1E2,
                            (cycle_3.contacts_mean-cycle_3.contacts_std),
                            (cycle_3.contacts_mean + cycle_3.contacts_std),
                            label=cycle_3.label, color='C5', alpha=1)

    lns_cycle_1 = ax_cycling.plot(cycle_1.strains * 1E2, cycle_1.contacts_mean, 'k-', zorder=3)
    ax_cycling.fill_between(cycle_1.strains * 1E2,
                            (cycle_1.contacts_mean-cycle_1.contacts_std),
                            (cycle_1.contacts_mean + cycle_1.contacts_std),
                            label=cycle_1.label, color='C4', alpha=1)

    lns_cycle_0 = ax_cycling.plot(SOC_0.strains * 1E2, SOC_0.contacts_mean, 'k-', zorder=3)
    ax_cycling.fill_between(SOC_0.strains * 1E2,
                        (SOC_0.contacts_mean-SOC_0.contacts_std),
                        (SOC_0.contacts_mean + SOC_0.contacts_std),
                        label=SOC_0.label, color='C1', alpha=1)

    ax_cycling.set_ylim(ymin=0)
    ax_cycling.set_xlim(xmin=-2.2, xmax=2.2)
    ax_cycling.set_ylabel('Particle contacts [-]')
    ax_cycling.set_xlabel('Strain [%]')
    ax_cycling.legend(loc='best', title='Charge cycles')

    fig_cycling.tight_layout()
    fname = fig_dir + 'cycling'
    plt.savefig(fname)

    print('Show plots')
    plt.show()

