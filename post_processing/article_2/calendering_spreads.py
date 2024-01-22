from force_model_impact_on_calendering.Bertil_calendering_pressure import local_data_gatherer, bertil_data_gatherer, \
    contact_counter_bertil
from force_model_impact_on_calendering.Bertil_calendering_pressure_multiple_simulations import \
    calendering_plot_processing
from force_model_impact_on_calendering.Bertil_mechanical_properties_multiple_runs import mech_plot_prop, stiffness_func

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import os
import matplotlib
import shutil
from scipy import interpolate

matplotlib.style.use('axel_style')


def calendering_spread_processing(sim_dir, sim_type, no_sims):
    simulation_dictionary = {}
    max_load_height = -1
    max_unload_height = -1
    min_calendering_height = 1e99
    for i in range(1, no_sims + 1):
        simulation_dictionary[i] = Calendering(sim_dir + str(i) + '/electrode_calendering_' + sim_type)
        temp_min_height = min(simulation_dictionary[i].calendering_surface_position)
        if min_calendering_height > temp_min_height:
            min_calendering_height = temp_min_height
        if simulation_dictionary[i].calendering_surface_position[0] > max_load_height:
            max_load_height = simulation_dictionary[i].calendering_surface_position[0]
        if simulation_dictionary[i].calendering_surface_position[-1] > max_unload_height:
            max_unload_height = simulation_dictionary[i].calendering_surface_position[-1]
    lin_space_size = 1000
    loading_array = np.linspace(min_calendering_height, max_load_height, lin_space_size)
    unloading_array = np.linspace(min_calendering_height, max_unload_height, lin_space_size)
    loading_pressure_matrix = np.zeros((no_sims, lin_space_size))
    unloading_pressure_matrix = np.zeros((no_sims, lin_space_size))
    for i in range(no_sims):
        loading_pressure_matrix[i, :] = simulation_dictionary[i + 1].loading_surface_pressure_interpol(loading_array)
        unloading_pressure_matrix[i, :] = simulation_dictionary[i + 1].unloading_surface_pressure_interpol(
            unloading_array)
    loading_surface_pressure_mean = np.mean(loading_pressure_matrix, axis=0)
    loading_surface_pressure_std = np.std(loading_pressure_matrix, axis=0)
    unloading_surface_pressure_mean = np.mean(unloading_pressure_matrix, axis=0)
    unloading_surface_pressure_std = np.std(unloading_pressure_matrix, axis=0)
    return loading_array, loading_surface_pressure_mean, loading_surface_pressure_std, \
           unloading_array, unloading_surface_pressure_mean, unloading_surface_pressure_std


class Calendering:
    def __init__(self, sim_dir):
        self.sim_dir = sim_dir

        self.time, self.calendering_surface_pressure, self.bottom_surface_pressure, self.calendering_surface_position, \
            self.sxx, self.syy, self.szz, self.tau_xy, self.tau_xz, self.tau_yx, self.tau_yz, self.tau_zx,\
            self.tau_zy, self.ke = calendering_plot_processing(sim_dir)
        i = 1
        while self.calendering_surface_position[i - 1] >= self.calendering_surface_position[i]:
            i += 1
        self.min_surface_height_index = i - 1
        self.loading_surface_position = self.calendering_surface_position[:self.min_surface_height_index]
        self.loading_surface_pressure = self.calendering_surface_pressure[:self.min_surface_height_index]
        self.loading_surface_pressure_interpol = interpolate.interp1d(self.loading_surface_position,
                                                                      self.loading_surface_pressure,
                                                                      bounds_error=False, fill_value=0)
        self.unloading_surface_position = self.calendering_surface_position[self.min_surface_height_index:]
        self.unloading_surface_pressure = self.calendering_surface_pressure[self.min_surface_height_index:]
        self.unloading_surface_pressure_interpol = interpolate.interp1d(self.unloading_surface_position,
                                                                        self.unloading_surface_pressure,
                                                                        bounds_error=False, fill_value=0)


class Simulation:
    def __init__(self, sim_dir, sim_type, no_sims, label='template'):
        self.sim_dir = sim_dir
        self.sim_type = sim_type
        self.no_sims = no_sims
        self.loading_surface_position, self.loading_surface_pressure_mean, self.loading_surface_pressure_std, \
            self.unloading_surface_position, self.unloading_surface_pressure_mean, self.unloading_surface_pressure_std \
            = calendering_spread_processing(sim_dir, sim_type, no_sims)
        self.label = label


if __name__ == '__main__':

    simulation_directory_SN_101 = '/scratch/users/axlun/DEMsim/results/article_2/final_runs_2/SN_101/'
    SN_101 = Simulation(simulation_directory_SN_101, 'hertz', 4, 'Reference')

    simulation_directory_SN_101P = '/scratch/users/axlun/DEMsim/results/article_2/final_runs_2/SN_101P/'
    SN_101P = Simulation(simulation_directory_SN_101P, 'hertz', 4, 'New PSD')

    simulation_directory_SN_111 = '/scratch/users/axlun/DEMsim/results/article_2/final_runs_2/SN_111/'
    SN_111 = Simulation(simulation_directory_SN_111, 'hertz', 4, 'SN_111')

    simulation_directory_SN_200 = '/scratch/users/axlun/DEMsim/results/article_2/final_runs_2/SN_200/'
    SN_200 = Simulation(simulation_directory_SN_200, 'el_pl_binder_el_pl_particle', 4, 'SN_200')

    simulation_directory_SN_201 = '/scratch/users/axlun/DEMsim/results/article_2/final_runs_2/SN_201/'
    SN_201 = Simulation(simulation_directory_SN_201, 'el_pl_binder_el_pl_particle', 4, 'Material 1')

    simulation_directory_SN_202 = '/scratch/users/axlun/DEMsim/results/article_2/final_runs_2/SN_202/'
    SN_202 = Simulation(simulation_directory_SN_202, 'el_pl_binder_el_pl_particle', 4, 'SN_202')

    # simulation_directory_SN_211 = '/scratch/users/axlun/DEMsim/results/article_2/final_runs_2/SN_211/'
    # SN_211 = Simulation(simulation_directory_SN_211, 'el_pl_binder_el_pl_particle', 4, 'SN_211')

    simulation_directory_SN_301 = '/scratch/users/axlun/DEMsim/results/article_2/final_runs_2/SN_301/'
    SN_301 = Simulation(simulation_directory_SN_301, 'el_pl_binder_el_pl_particle', 4, 'Material 2')

    # ==FIGURE SAVE DIR=================================================================================================
    fig_dir = 'C:/temp/figures/calendering_spreads/'
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

    # =FIGURE PARTICLE SIZE=============================================================================================
    fig_particle_size, ax_particle_size = plt.subplots()
    ax_particle_size.plot(SN_101.loading_surface_position * 1E2,
                          SN_101.loading_surface_pressure_mean * 1E-6,
                          'k-', zorder=1)
    ax_particle_size.plot(SN_101.unloading_surface_position * 1E2,
                          SN_101.unloading_surface_pressure_mean * 1E-6,
                          'k-', zorder=1)
    ax_particle_size.fill_between(
        SN_101.loading_surface_position * 1E2,
        (SN_101.loading_surface_pressure_mean - SN_101.loading_surface_pressure_std) * 1E-6,
        (SN_101.loading_surface_pressure_mean + SN_101.loading_surface_pressure_std) * 1E-6,
        color='C0', label=SN_101.label, alpha=1)
    ax_particle_size.fill_between(
        SN_101.unloading_surface_position * 1E2,
        (SN_101.unloading_surface_pressure_mean - SN_101.unloading_surface_pressure_std) * 1E-6,
        (SN_101.unloading_surface_pressure_mean + SN_101.unloading_surface_pressure_std) * 1E-6,
        color='C0', alpha=1)

    ax_particle_size.plot(SN_101P.loading_surface_position * 1E2,
                          SN_101P.loading_surface_pressure_mean * 1E-6,
                          'k-', zorder=2)
    ax_particle_size.plot(SN_101P.unloading_surface_position * 1E2,
                          SN_101P.unloading_surface_pressure_mean * 1E-6,
                          'k-', zorder=2)
    ax_particle_size.fill_between(
        SN_101P.loading_surface_position * 1E2,
        (SN_101P.loading_surface_pressure_mean - SN_101P.loading_surface_pressure_std) * 1E-6,
        (SN_101P.loading_surface_pressure_mean + SN_101P.loading_surface_pressure_std) * 1E-6,
        color='C1', zorder=2, label=SN_101P.label, alpha=1)
    ax_particle_size.fill_between(
        SN_101P.unloading_surface_position * 1E2,
        (SN_101P.unloading_surface_pressure_mean - SN_101P.unloading_surface_pressure_std) * 1E-6,
        (SN_101P.unloading_surface_pressure_mean + SN_101P.unloading_surface_pressure_std) * 1E-6,
        color='C1', zorder=2, alpha=1)

    ax_particle_size.set_xlim(xmin=104.8, xmax=120)
    ax_particle_size.xaxis.set_major_locator(MultipleLocator(5))
    ax_particle_size.set_ylim(ymin=0)
    ax_particle_size.set_ylabel('Calendering surface pressure [MPa]')
    ax_particle_size.set_xlabel('Calendering surface height [µm]')
    fig_particle_size.tight_layout()
    ax_particle_size.legend(loc='best')
    fname = fig_dir + 'size_distribution'
    plt.savefig(fname)
    # =FIGURE CONTACT MODEL=============================================================================================
    fig_calendering, ax_calendering = plt.subplots()
    """
    ax_calendering.plot(SN_111.loading_surface_position * 1E2,
                        SN_111.loading_surface_pressure_mean * 1E-6,
                        'k-', zorder=1)
    ax_calendering.plot(SN_111.unloading_surface_position * 1E2,
                        SN_111.unloading_surface_pressure_mean * 1E-6,
                        'k-', zorder=1)
    ax_calendering.fill_between(SN_111.loading_surface_position * 1E2,
                                (SN_111.loading_surface_pressure_mean - SN_111.loading_surface_pressure_std)*1E-6,
                                (SN_111.loading_surface_pressure_mean + SN_111.loading_surface_pressure_std)*1E-6,
                                color='C1', label=SN_111.label, alpha=1)
    ax_calendering.fill_between(SN_111.unloading_surface_position*1E2,
                                (SN_111.unloading_surface_pressure_mean - SN_111.unloading_surface_pressure_std) * 1E-6,
                                (SN_111.unloading_surface_pressure_mean + SN_111.unloading_surface_pressure_std) * 1E-6,
                                color='C1', alpha=1)
    """
    """"
    ax_calendering.plot(SN_211.loading_surface_position * 1E2,
                        SN_211.loading_surface_pressure_mean * 1E-6,
                        'k-', zorder=1)
    ax_calendering.plot(SN_211.unloading_surface_position * 1E2,
                        SN_211.unloading_surface_pressure_mean * 1E-6,
                        'k-', zorder=1)
    ax_calendering.fill_between(SN_211.loading_surface_position * 1E2,
                                (SN_211.loading_surface_pressure_mean - SN_211.loading_surface_pressure_std)*1E-6,
                                (SN_211.loading_surface_pressure_mean + SN_211.loading_surface_pressure_std)*1E-6,
                                color='C2', label=SN_211.label, alpha=1)
    ax_calendering.fill_between(SN_211.unloading_surface_position*1E2,
                                (SN_211.unloading_surface_pressure_mean - SN_211.unloading_surface_pressure_std) * 1E-6,
                                (SN_211.unloading_surface_pressure_mean + SN_211.unloading_surface_pressure_std) * 1E-6,
                                color='C2', alpha=1)
    """
    ax_calendering.plot(SN_101.loading_surface_position * 1E2,
                        SN_101.loading_surface_pressure_mean * 1E-6,
                        'k-', zorder=3)
    ax_calendering.plot(SN_101.unloading_surface_position * 1E2,
                        SN_101.unloading_surface_pressure_mean * 1E-6,
                        'k-', zorder=3)
    ax_calendering.fill_between(SN_101.loading_surface_position * 1E2,
                                (SN_101.loading_surface_pressure_mean - SN_101.loading_surface_pressure_std) * 1E-6,
                                (SN_101.loading_surface_pressure_mean + SN_101.loading_surface_pressure_std) * 1E-6,
                                color='C0', zorder=3, label=SN_101.label, alpha=1)
    ax_calendering.fill_between(SN_101.unloading_surface_position * 1E2,
                                (SN_101.unloading_surface_pressure_mean - SN_101.unloading_surface_pressure_std) * 1E-6,
                                (SN_101.unloading_surface_pressure_mean + SN_101.unloading_surface_pressure_std) * 1E-6,
                                color='C0', zorder=3, alpha=1)

    ax_calendering.plot(SN_201.loading_surface_position * 1E2,
                        SN_201.loading_surface_pressure_mean * 1E-6,
                        'k-', zorder=2)
    ax_calendering.plot(SN_201.unloading_surface_position * 1E2,
                        SN_201.unloading_surface_pressure_mean * 1E-6,
                        'k-', zorder=2)
    ax_calendering.fill_between(SN_201.loading_surface_position * 1E2,
                                (SN_201.loading_surface_pressure_mean - SN_201.loading_surface_pressure_std) * 1E-6,
                                (SN_201.loading_surface_pressure_mean + SN_201.loading_surface_pressure_std) * 1E-6,
                                color='C2', zorder=2, label=SN_201.label, alpha=1)
    ax_calendering.fill_between(SN_201.unloading_surface_position * 1E2,
                                (SN_201.unloading_surface_pressure_mean - SN_201.unloading_surface_pressure_std) * 1E-6,
                                (SN_201.unloading_surface_pressure_mean + SN_201.unloading_surface_pressure_std) * 1E-6,
                                color='C2', zorder=2, alpha=1)

    ax_calendering.plot(SN_301.loading_surface_position * 1E2,
                        SN_301.loading_surface_pressure_mean * 1E-6,
                        'k-', zorder=1)
    ax_calendering.plot(SN_301.unloading_surface_position * 1E2,
                        SN_301.unloading_surface_pressure_mean * 1E-6,
                        'k-', zorder=1)
    ax_calendering.fill_between(SN_301.loading_surface_position * 1E2,
                                (SN_301.loading_surface_pressure_mean - SN_301.loading_surface_pressure_std) * 1E-6,
                                (SN_301.loading_surface_pressure_mean + SN_301.loading_surface_pressure_std) * 1E-6,
                                color='C3', zorder=1, label=SN_301.label, alpha=1)
    ax_calendering.fill_between(SN_301.unloading_surface_position * 1E2,
                                (SN_301.unloading_surface_pressure_mean - SN_301.unloading_surface_pressure_std) * 1E-6,
                                (SN_301.unloading_surface_pressure_mean + SN_301.unloading_surface_pressure_std) * 1E-6,
                                color='C3', zorder=1, alpha=1)

    ax_calendering.set_xlim(xmin=104.8, xmax=120)
    ax_calendering.xaxis.set_major_locator(MultipleLocator(5))
    ax_calendering.set_ylim(ymin=0)
    ax_calendering.set_ylabel('Calendering surface pressure [MPa]')
    ax_calendering.set_xlabel('Calendering surface height [µm]')
    fig_calendering.tight_layout()

    # handles_sim = lns_fill_SN_101 + lns_fill_SN_201 + lns_fill_SN_301
    # labels_sim = [l.get_label() for l in handles_sim]
    # plt.legend(handles_sim, labels_sim, loc='upper right')#, title='Simulations')
    ax_calendering.legend(loc='best')
    fname = fig_dir + 'contact_model'
    plt.savefig(fname)

    # =FIGURE CALENDERING DEGREE======================================================================================
    fig_calendering_degree, ax_calendering_degree = plt.subplots()
    ax_calendering_degree.plot(SN_202.loading_surface_position*1E2,
                               SN_202.loading_surface_pressure_mean*1E-6,
                               'k-', zorder=1)
    ax_calendering_degree.plot(SN_202.unloading_surface_position*1E2,
                               SN_202.unloading_surface_pressure_mean*1E-6,
                               'k-', zorder=1)
    ax_calendering_degree.fill_between(SN_202.loading_surface_position*1E2,
                               (SN_202.loading_surface_pressure_mean-SN_202.loading_surface_pressure_std)*1E-6,
                               (SN_202.loading_surface_pressure_mean+SN_202.loading_surface_pressure_std)*1E-6,
                               color='C2', label=SN_202.label, alpha=1)
    ax_calendering_degree.fill_between(SN_202.unloading_surface_position*1E2,
                               (SN_202.unloading_surface_pressure_mean-SN_202.unloading_surface_pressure_std)*1E-6,
                               (SN_202.unloading_surface_pressure_mean+SN_202.unloading_surface_pressure_std)*1E-6,
                               color='C2', alpha=1)

    ax_calendering_degree.plot(SN_201.loading_surface_position * 1E2,
                               SN_201.loading_surface_pressure_mean * 1E-6,
                               'k-', zorder=1)
    ax_calendering_degree.plot(SN_201.unloading_surface_position * 1E2,
                               SN_201.unloading_surface_pressure_mean * 1E-6,
                               'k-', zorder=1)
    ax_calendering_degree.fill_between(
        SN_201.loading_surface_position * 1E2,
        (SN_201.loading_surface_pressure_mean - SN_201.loading_surface_pressure_std) * 1E-6,
        (SN_201.loading_surface_pressure_mean + SN_201.loading_surface_pressure_std) * 1E-6,
        color='C1', label=SN_201.label, alpha=1)
    ax_calendering_degree.fill_between(
        SN_201.unloading_surface_position * 1E2,
        (SN_201.unloading_surface_pressure_mean - SN_201.unloading_surface_pressure_std) * 1E-6,
        (SN_201.unloading_surface_pressure_mean + SN_201.unloading_surface_pressure_std) * 1E-6,
        color='C1', alpha=1)

    ax_calendering_degree.plot(SN_200.loading_surface_position * 1E2,
                               SN_200.loading_surface_pressure_mean * 1E-6,
                               'k-', zorder=1)
    ax_calendering_degree.plot(SN_200.unloading_surface_position * 1E2,
                               SN_200.unloading_surface_pressure_mean * 1E-6,
                               'k-', zorder=1)
    ax_calendering_degree.fill_between(
        SN_200.loading_surface_position * 1E2,
        (SN_200.loading_surface_pressure_mean - SN_200.loading_surface_pressure_std) * 1E-6,
        (SN_200.loading_surface_pressure_mean + SN_200.loading_surface_pressure_std) * 1E-6,
        color='C0', label=SN_200.label, alpha=1)
    ax_calendering_degree.fill_between(
        SN_200.unloading_surface_position * 1E2,
        (SN_200.unloading_surface_pressure_mean - SN_200.unloading_surface_pressure_std) * 1E-6,
        (SN_200.unloading_surface_pressure_mean + SN_200.unloading_surface_pressure_std) * 1E-6,
        color='C0', alpha=1)

    ax_calendering_degree.set_xlim(xmin=102.2, xmax=120)
    ax_calendering_degree.xaxis.set_major_locator(MultipleLocator(5))
    ax_calendering_degree.set_ylim(ymin=0)
    ax_calendering_degree.set_ylabel('Calendering surface pressure [MPa]')
    ax_calendering_degree.set_xlabel('Calendering surface height [µm]')
    fig_calendering_degree.tight_layout()
    ax_calendering_degree.legend(loc='best')
    fname = fig_dir + 'calendering_degree'
    plt.savefig(fname)
    plt.show()
