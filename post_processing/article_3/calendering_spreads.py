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


def calendering_spread_processing(sim_dir, no_sims):
    simulation_dictionary = {}
    max_load_height = -1
    max_unload_height = -1
    min_calendering_height = 1e99
    for i in range(1, no_sims + 1):
        simulation_dictionary[i] = Calendering(sim_dir + str(i) + '/swelling_electrode_calendering')
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
        self.min_surface_height_index = i
        self.loading_surface_position = self.calendering_surface_position[:self.min_surface_height_index]
        self.loading_surface_pressure = self.calendering_surface_pressure[:self.min_surface_height_index]
        self.loading_surface_pressure_interpol = interpolate.interp1d(self.loading_surface_position,
                                                                      self.loading_surface_pressure,
                                                                      bounds_error=False, fill_value=0)
        self.unloading_surface_position = self.calendering_surface_position[(self.min_surface_height_index-1):]
        self.unloading_surface_pressure = self.calendering_surface_pressure[(self.min_surface_height_index-1):]
        self.unloading_surface_pressure_interpol = interpolate.interp1d(self.unloading_surface_position,
                                                                        self.unloading_surface_pressure,
                                                                        bounds_error=False, fill_value=0)


class Simulation:
    def __init__(self, sim_dir, no_sims, label='template'):
        self.sim_dir = sim_dir
        self.no_sims = no_sims
        self.loading_surface_position, self.loading_surface_pressure_mean, self.loading_surface_pressure_std, \
            self.unloading_surface_position, self.unloading_surface_pressure_mean, self.unloading_surface_pressure_std \
            = calendering_spread_processing(sim_dir, no_sims)
        self.label = label


if __name__ == '__main__':

    calendering = Simulation('/scratch/users/axlun/DEMsim/results/article_3/final_runs/', 4, 'Reference')

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
    ax_particle_size.plot(calendering.loading_surface_position * 1E2,
                          calendering.loading_surface_pressure_mean * 1E-6,
                          'k-', zorder=1)
    ax_particle_size.plot(calendering.unloading_surface_position * 1E2,
                          calendering.unloading_surface_pressure_mean * 1E-6,
                          'k-', zorder=1)
    ax_particle_size.fill_between(
        calendering.loading_surface_position * 1E2,
        (calendering.loading_surface_pressure_mean - calendering.loading_surface_pressure_std) * 1E-6,
        (calendering.loading_surface_pressure_mean + calendering.loading_surface_pressure_std) * 1E-6,
        color='C0', label=calendering.label, alpha=1)
    ax_particle_size.fill_between(
        calendering.unloading_surface_position * 1E2,
        (calendering.unloading_surface_pressure_mean - calendering.unloading_surface_pressure_std) * 1E-6,
        (calendering.unloading_surface_pressure_mean + calendering.unloading_surface_pressure_std) * 1E-6,
        color='C0', alpha=1)

    ax_particle_size.set_xlim(xmin=104.8, xmax=120)
    ax_particle_size.xaxis.set_major_locator(MultipleLocator(5))
    ax_particle_size.set_ylim(ymin=0)
    ax_particle_size.set_ylabel('Calendering surface pressure [MPa]')
    ax_particle_size.set_xlabel('Calendering surface height [µm]')
    fig_particle_size.tight_layout()
    ax_particle_size.legend(loc='best')
    fname = fig_dir + 'size_distribution'
    plt.savefig(fname)
    plt.show()
