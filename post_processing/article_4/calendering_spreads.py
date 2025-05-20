#
# Created by Axel on 2025-05-08
# Adapted from article_3/calendering_spreads.py
#
from force_model_impact_on_calendering.Bertil_calendering_pressure import local_data_gatherer, bertil_data_gatherer, \
    contact_counter_bertil
from force_model_impact_on_calendering.Bertil_calendering_pressure_multiple_simulations import \
    calendering_plot_processing
from force_model_impact_on_calendering.Bertil_fractured_particles import fractured_particle_gatherer
from force_model_impact_on_calendering.Bertil_fractured_binder_contacts import fractured_binder_counter
from force_model_impact_on_calendering.Bertil_mechanical_properties_multiple_runs import mech_plot_prop, stiffness_func

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import os
import matplotlib
import shutil
from scipy import interpolate

matplotlib.style.use('axel_style')


# Removes zero values in the beginning and end of the matrix rows and replaces them with the first or last non-zero
# value of the row
def fill_matrix_edges(matrix):
    result = matrix.copy()
    rows, cols = result.shape
    for i in range(rows):
        row = result[i]
        nz = np.flatnonzero(row)
        if nz.size == 0:
            continue
        row[:nz[0]] = row[nz[0]]
        row[nz[-1] + 1:] = row[nz[-1]]
    return result


def calendering_spread_processing(sim_case, no_sims):
    simulation_dictionary = {}
    max_load_height = -1
    max_unload_height = -1
    min_calendering_height = 1e99
    for i in range(1, no_sims + 1):
        print('/scratch/users/axlun/DEMsim/results/article_4/final_runs/' + str(i) + '/' + sim_case + '/')
        simulation_dictionary[i] = Calendering('/scratch/users/axlun/DEMsim/results/article_4/final_runs/'
                                               + str(i) + '/' + sim_case + '/electrode_calendering')
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
    loading_binder_contacts_matrix = np.zeros((no_sims, lin_space_size))
    loading_initialised_binder_contacts_matrix = np.zeros((no_sims, lin_space_size))
    loading_fractured_binders_matrix = np.zeros((no_sims, lin_space_size))
    loading_particle_contacts_matrix = np.zeros((no_sims, lin_space_size))
    loading_fractured_particles_matrix = np.zeros((no_sims, lin_space_size))
    unloading_pressure_matrix = np.zeros((no_sims, lin_space_size))
    unloading_binder_contacts_matrix = np.zeros((no_sims, lin_space_size))
    unloading_initialised_binder_contacts_matrix = np.zeros((no_sims, lin_space_size))
    unloading_fractured_binders_matrix = np.zeros((no_sims, lin_space_size))
    unloading_particle_contacts_matrix = np.zeros((no_sims, lin_space_size))
    unloading_fractured_particles_matrix = np.zeros((no_sims, lin_space_size))

    for i in range(no_sims):
        loading_pressure_matrix[i, :] = simulation_dictionary[i + 1].loading_surface_pressure_interpol(loading_array)
        loading_initialised_binder_contacts_matrix[i, :] = \
            simulation_dictionary[i + 1].loading_initialised_binder_contacts_interpol(loading_array)
        loading_binder_contacts_matrix[i, :] = \
            simulation_dictionary[i + 1].loading_active_binder_contacts_interpol(loading_array)
        loading_fractured_binders_matrix[i, :] = \
            simulation_dictionary[i + 1].loading_fractured_binder_interpol(loading_array)
        loading_particle_contacts_matrix[i, :] = \
            simulation_dictionary[i + 1].loading_particle_contacts_interpol(loading_array)
        loading_fractured_particles_matrix[i, :] = \
            simulation_dictionary[i + 1].loading_fractured_particles_interpol(loading_array)

        unloading_pressure_matrix[i, :] = \
            simulation_dictionary[i + 1].unloading_surface_pressure_interpol(unloading_array)
        unloading_initialised_binder_contacts_matrix[i, :] = \
            simulation_dictionary[i + 1].unloading_initialised_binder_contacts_interpol(unloading_array)
        unloading_binder_contacts_matrix[i, :] = \
            simulation_dictionary[i + 1].unloading_active_binder_contacts_interpol(unloading_array)
        unloading_fractured_binders_matrix[i, :] = \
            simulation_dictionary[i + 1].unloading_fracture_binder_interpol(unloading_array)
        unloading_particle_contacts_matrix[i, :] =\
            simulation_dictionary[i + 1].unloading_particle_contacts_interpol(unloading_array)
        unloading_fractured_particles_matrix[i, :] = \
            simulation_dictionary[i + 1].unloading_fracture_particles_interpol(unloading_array)

    # =REMOVE NON ZERO VALUES IN BEGINNING AND END OF MATRIX============================================================
    loading_binder_contacts_matrix = fill_matrix_edges(loading_binder_contacts_matrix)
    loading_initialised_binder_contacts_matrix = fill_matrix_edges(loading_initialised_binder_contacts_matrix)
    loading_fractured_binders_matrix = fill_matrix_edges(loading_fractured_binders_matrix)
    loading_particle_contacts_matrix = fill_matrix_edges(loading_particle_contacts_matrix)
    loading_fractured_particles_matrix = fill_matrix_edges(loading_fractured_particles_matrix)

    unloading_binder_contacts_matrix = fill_matrix_edges(unloading_binder_contacts_matrix)
    unloading_initialised_binder_contacts_matrix = fill_matrix_edges(unloading_initialised_binder_contacts_matrix)
    unloading_fractured_binders_matrix = fill_matrix_edges(unloading_fractured_binders_matrix)
    unloading_particle_contacts_matrix = fill_matrix_edges(unloading_particle_contacts_matrix)
    unloading_fractured_particles_matrix = fill_matrix_edges(unloading_fractured_particles_matrix)

    # =NORMALISATION OF FRACTURE MATRIX=================================================================================
    loading_fractured_binders_matrix = loading_fractured_binders_matrix / loading_initialised_binder_contacts_matrix
    unloading_fractured_binders_matrix = unloading_fractured_binders_matrix/unloading_initialised_binder_contacts_matrix
    # loading_fractured_binders_matrix = loading_fractured_binders_matrix / loading_binder_contacts_matrix
    # unloading_fractured_binders_matrix = unloading_fractured_binders_matrix / unloading_binder_contacts_matrix

    loading_fractured_particles_matrix = loading_fractured_particles_matrix / 5000
    unloading_fractured_particles_matrix = unloading_fractured_particles_matrix / 5000
    # ==================================================================================================================

    loading_surface_pressure_mean = np.mean(loading_pressure_matrix, axis=0)
    loading_surface_pressure_std = np.std(loading_pressure_matrix, axis=0)

    loading_binder_contacts_mean = np.mean(loading_binder_contacts_matrix, axis=0)
    loading_binder_contacts_std = np.std(loading_binder_contacts_matrix, axis=0)

    loading_fractured_binder_mean = np.mean(loading_fractured_binders_matrix, axis=0)
    loading_fractured_binder_std = np.std(loading_fractured_binders_matrix, axis=0)

    loading_particle_contacts_mean = np.mean(loading_particle_contacts_matrix, axis=0)
    loading_particle_contacts_std = np.std(loading_particle_contacts_matrix, axis=0)

    loading_fractured_particles_mean = np.mean(loading_fractured_particles_matrix, axis=0)
    loading_fractured_particles_std = np.std(loading_fractured_particles_matrix, axis=0)

    unloading_surface_pressure_mean = np.mean(unloading_pressure_matrix, axis=0)
    unloading_surface_pressure_std = np.std(unloading_pressure_matrix, axis=0)

    unloading_binder_contacts_mean = np.mean(unloading_binder_contacts_matrix, axis=0)
    unloading_binder_contacts_std = np.std(unloading_binder_contacts_matrix, axis=0)

    unloading_fractured_binder_mean = np.mean(unloading_fractured_binders_matrix, axis=0)
    unloading_fractured_binder_std = np.std(unloading_fractured_binders_matrix, axis=0)

    unloading_particle_contacts_mean = np.mean(unloading_particle_contacts_matrix, axis=0)
    unloading_particle_contacts_std = np.std(unloading_particle_contacts_matrix, axis=0)

    unloading_fractured_particles_mean = np.mean(unloading_fractured_particles_matrix, axis=0)
    unloading_fractured_particles_std = np.std(unloading_fractured_particles_matrix, axis=0)

    return loading_array, loading_surface_pressure_mean, loading_surface_pressure_std, loading_pressure_matrix, \
           loading_particle_contacts_mean, loading_particle_contacts_std, loading_particle_contacts_matrix, \
           loading_binder_contacts_mean, loading_binder_contacts_std, loading_binder_contacts_matrix, \
           loading_fractured_particles_mean, loading_fractured_particles_std, loading_fractured_particles_matrix, \
           loading_fractured_binder_mean, loading_fractured_binder_std, loading_fractured_binders_matrix, \
           unloading_array, unloading_surface_pressure_mean, unloading_surface_pressure_std, unloading_pressure_matrix, \
           unloading_particle_contacts_mean, unloading_particle_contacts_std, unloading_particle_contacts_matrix, \
           unloading_binder_contacts_mean, unloading_binder_contacts_std, unloading_binder_contacts_matrix, \
           unloading_fractured_particles_mean, unloading_fractured_particles_std, unloading_fractured_particles_matrix, \
           unloading_fractured_binder_mean, unloading_fractured_binder_std, unloading_fractured_binders_matrix


class Calendering:
    def __init__(self, sim_dir):
        self.sim_dir = sim_dir

        self.time, self.calendering_surface_pressure, self.bottom_surface_pressure, self.calendering_surface_position, \
        self.sxx, self.syy, self.szz, self.tau_xy, self.tau_xz, self.tau_yx, self.tau_yz, self.tau_zx, self.tau_zy, \
        self.ke = calendering_plot_processing(sim_dir)

        # Also get the number of contacts here
        self.contacts_time, self.particle_contact_vec, self.binder_contact_vec, self.binder_particle_contact_vec \
            = contact_counter_bertil(sim_dir)

        # As well as the number of fractures
        self.frac_binder_time, self.initialised_binder_contact_vec, self.active_binder_contact_vec, \
        self.fractured_binder_contact_vec = fractured_binder_counter(sim_dir)

        self.frac_particle_time, self.fractured_particles, self.fracture_particle_array = \
            fractured_particle_gatherer(sim_dir)

        i = 1
        while self.calendering_surface_position[i - 1] >= self.calendering_surface_position[i]:
            i += 1
        self.min_surface_height_index = i
        self.loading_surface_position = self.calendering_surface_position[:self.min_surface_height_index]
        self.loading_surface_pressure = self.calendering_surface_pressure[:self.min_surface_height_index]
        self.loading_surface_pressure_interpol = interpolate.interp1d(self.loading_surface_position,
                                                                      self.loading_surface_pressure,
                                                                      bounds_error=False, fill_value=0)

        # Split and interpolate contact vectors and fracture vectors for loading here===================================
        self.loading_initialised_binder_contacts = self.initialised_binder_contact_vec[:self.min_surface_height_index]
        self.loading_initialised_binder_contacts_interpol = interpolate.interp1d(self.loading_surface_position,
                                                                                 self.loading_initialised_binder_contacts,
                                                                                 bounds_error=False, fill_value=0)

        self.loading_active_binder_contacts = self.active_binder_contact_vec[:self.min_surface_height_index]
        self.loading_active_binder_contacts_interpol = interpolate.interp1d(self.loading_surface_position,
                                                                            self.loading_active_binder_contacts,
                                                                            bounds_error=False, fill_value=0)

        self.loading_fractured_binder = self.fractured_binder_contact_vec[:self.min_surface_height_index]
        self.loading_fractured_binder_interpol = interpolate.interp1d(self.loading_surface_position,
                                                                      self.loading_fractured_binder,
                                                                      bounds_error=False, fill_value=0)

        self.loading_particle_contacts = self.particle_contact_vec[:self.min_surface_height_index]
        self.loading_particle_contacts_interpol = interpolate.interp1d(self.loading_surface_position,
                                                                       self.loading_particle_contacts,
                                                                       bounds_error=False, fill_value=0)

        self.loading_fractured_particles = self.fractured_particles[:self.min_surface_height_index]
        self.loading_fractured_particles_interpol = interpolate.interp1d(self.loading_surface_position,
                                                                         self.loading_fractured_particles,
                                                                         bounds_error=False, fill_value=0)
        # ==============================================================================================================

        self.unloading_surface_position = self.calendering_surface_position[(self.min_surface_height_index - 1):]
        self.unloading_surface_pressure = self.calendering_surface_pressure[(self.min_surface_height_index - 1):]
        self.unloading_surface_pressure_interpol = interpolate.interp1d(self.unloading_surface_position,
                                                                        self.unloading_surface_pressure,
                                                                        bounds_error=False, fill_value=0)

        # Split and interpolate contact vectors and fracture vectors for unloading here=================================
        self.unloading_initialised_binder_contacts = self.initialised_binder_contact_vec[
                                                     (self.min_surface_height_index - 1):]
        self.unloading_initialised_binder_contacts_interpol = interpolate.interp1d(self.unloading_surface_position,
                                                                                   self.unloading_initialised_binder_contacts,
                                                                                   bounds_error=False, fill_value=0)

        self.unloading_active_binder_contacts = self.active_binder_contact_vec[(self.min_surface_height_index - 1):]
        self.unloading_active_binder_contacts_interpol = interpolate.interp1d(self.unloading_surface_position,
                                                                              self.unloading_active_binder_contacts,
                                                                              bounds_error=False, fill_value=0)

        self.unloading_fracture_binder = self.fractured_binder_contact_vec[(self.min_surface_height_index - 1):]
        self.unloading_fracture_binder_interpol = interpolate.interp1d(self.unloading_surface_position,
                                                                       self.unloading_fracture_binder,
                                                                       bounds_error=False, fill_value=0)

        self.unloading_particle_contacts = self.particle_contact_vec[(self.min_surface_height_index - 1):]
        self.unloading_particle_contacts_interpol = interpolate.interp1d(self.unloading_surface_position,
                                                                         self.unloading_particle_contacts,
                                                                         bounds_error=False, fill_value=0)

        self.unloading_fracture_particles = self.fractured_particles[(self.min_surface_height_index - 1):]
        self.unloading_fracture_particles_interpol = interpolate.interp1d(self.unloading_surface_position,
                                                                          self.unloading_fracture_particles,
                                                                          bounds_error=False, fill_value=0)
        # ===============================================================================================================


class Simulation:
    def __init__(self, sim_case, no_sims, label='template'):
        self.sim_dir = sim_case
        self.no_sims = no_sims
        self.loading_surface_position, self.loading_surface_pressure_mean, self.loading_surface_pressure_std, \
        self.loading_surface_pressure_matrix, self.loading_particle_contacts_mean, self.loading_particle_contacts_std, \
        self.loading_particle_contacts_matrix, self.loading_binder_contacts_mean, self.loading_binder_contacts_std, \
        self.loading_binder_contacts_matrix, self.loading_fractured_particles_mean, \
        self.loading_fractured_particles_std, self.loading_fractured_particles_matrix, \
        self.loading_fractured_binders_mean, self.loading_fractured_binders_std, self.loading_fractured_binders_matrix, \
        self.unloading_surface_position, self.unloading_surface_pressure_mean, self.unloading_surface_pressure_std, \
        self.unloading_surface_pressure_matrix, self.unloading_particle_contacts_mean, \
        self.unloading_particle_contacts_std, self.unloading_particle_contacts_matrix, \
        self.unloading_binder_contacts_mean, self.unloading_binder_contacts_std, \
        self.unloading_binder_contacts_matrix, self.unloading_fractured_particles_mean, \
        self.unloading_fractured_particles_std, self.unloading_fractured_particles_matrix, \
        self.unloading_fractured_binders_mean, self.unloading_fractured_binders_std, \
        self.unloading_fractured_binders_matrix = calendering_spread_processing(sim_case, no_sims)
        self.label = label


if __name__ == '__main__':
    # ==================================================================================================================
    # TODO: * How are binder contacts counted?
    #       * How are particle contacts counted?
    #           - Are combinded contacts captured?
    # ==================================================================================================================

    # ==FIGURE SAVE DIR=================================================================================================
    fig_dir = 'C:/temp/figures/article_4/final_runs/electrode_calendering_spreads/'
    try:
        shutil.rmtree(fig_dir)
        os.makedirs(fig_dir)
    except FileNotFoundError:
        print('Directory does not exist')
        os.makedirs(fig_dir)
    except:
        print('Directory could not be removed')
        quit()
    # ==================================================================================================================
    no_sims = 4
    ref = Simulation('ref', no_sims, 'Reference')
    sf_2e8 = Simulation('sf_2e8', no_sims, '$\sigma_{f}=200$ MPa')
    sf_1e8 = Simulation('sf_1e8', no_sims, '$\sigma_{f}=100$ MPa')
    sf_5e7 = Simulation('sf_5e7', no_sims, '$\sigma_{f}=50$ MPa')
    ########################################PARTICLE FRACTURE###########################################################
    # =FIGURE PARTICLE PRESSURE=========================================================================================
    fig_particle_pressure, ax_particle_pressure = plt.subplots()
    ax_particle_pressure.plot(ref.loading_surface_position * 1E2,
                              ref.loading_surface_pressure_mean * 1E-6,
                              'k-', zorder=1)
    ax_particle_pressure.plot(ref.unloading_surface_position * 1E2,
                              ref.unloading_surface_pressure_mean * 1E-6,
                              'k-', zorder=1)
    ax_particle_pressure.fill_between(
        ref.loading_surface_position * 1E2,
        (ref.loading_surface_pressure_mean - ref.loading_surface_pressure_std) * 1E-6,
        (ref.loading_surface_pressure_mean + ref.loading_surface_pressure_std) * 1E-6,
        color='C0', label=ref.label, alpha=1, zorder=1)
    ax_particle_pressure.fill_between(
        ref.unloading_surface_position * 1E2,
        (ref.unloading_surface_pressure_mean - ref.unloading_surface_pressure_std) * 1E-6,
        (ref.unloading_surface_pressure_mean + ref.unloading_surface_pressure_std) * 1E-6,
        color='C0', alpha=1, zorder=1)

    ax_particle_pressure.plot(sf_2e8.loading_surface_position * 1E2,
                              sf_2e8.loading_surface_pressure_mean * 1E-6,
                              'k-', zorder=2)
    ax_particle_pressure.plot(sf_2e8.unloading_surface_position * 1E2,
                              sf_2e8.unloading_surface_pressure_mean * 1E-6,
                              'k-', zorder=2)
    ax_particle_pressure.fill_between(
        sf_2e8.loading_surface_position * 1E2,
        (sf_2e8.loading_surface_pressure_mean - sf_2e8.loading_surface_pressure_std) * 1E-6,
        (sf_2e8.loading_surface_pressure_mean + sf_2e8.loading_surface_pressure_std) * 1E-6,
        color='C1', label=sf_2e8.label, alpha=1, zorder=2)
    ax_particle_pressure.fill_between(
        sf_2e8.unloading_surface_position * 1E2,
        (sf_2e8.unloading_surface_pressure_mean - sf_2e8.unloading_surface_pressure_std) * 1E-6,
        (sf_2e8.unloading_surface_pressure_mean + sf_2e8.unloading_surface_pressure_std) * 1E-6,
        color='C1', alpha=1, zorder=2)

    ax_particle_pressure.plot(sf_1e8.loading_surface_position * 1E2,
                              sf_1e8.loading_surface_pressure_mean * 1E-6,
                              'k-', zorder=3)
    ax_particle_pressure.plot(sf_1e8.unloading_surface_position * 1E2,
                              sf_1e8.unloading_surface_pressure_mean * 1E-6,
                              'k-', zorder=3)
    ax_particle_pressure.fill_between(
        sf_1e8.loading_surface_position * 1E2,
        (sf_1e8.loading_surface_pressure_mean - sf_1e8.loading_surface_pressure_std) * 1E-6,
        (sf_1e8.loading_surface_pressure_mean + sf_1e8.loading_surface_pressure_std) * 1E-6,
        color='C2', label=sf_1e8.label, alpha=1, zorder=3)
    ax_particle_pressure.fill_between(
        sf_1e8.unloading_surface_position * 1E2,
        (sf_1e8.unloading_surface_pressure_mean - sf_1e8.unloading_surface_pressure_std) * 1E-6,
        (sf_1e8.unloading_surface_pressure_mean + sf_1e8.unloading_surface_pressure_std) * 1E-6,
        color='C2', alpha=1, zorder=3)

    ax_particle_pressure.plot(sf_5e7.loading_surface_position * 1E2,
                              sf_5e7.loading_surface_pressure_mean * 1E-6,
                              'k-', zorder=4)
    ax_particle_pressure.plot(sf_5e7.unloading_surface_position * 1E2,
                              sf_5e7.unloading_surface_pressure_mean * 1E-6,
                              'k-', zorder=4)
    ax_particle_pressure.fill_between(
        sf_5e7.loading_surface_position * 1E2,
        (sf_5e7.loading_surface_pressure_mean - sf_5e7.loading_surface_pressure_std) * 1E-6,
        (sf_5e7.loading_surface_pressure_mean + sf_5e7.loading_surface_pressure_std) * 1E-6,
        color='C3', label=sf_5e7.label, alpha=1, zorder=4)
    ax_particle_pressure.fill_between(
        sf_5e7.unloading_surface_position * 1E2,
        (sf_5e7.unloading_surface_pressure_mean - sf_5e7.unloading_surface_pressure_std) * 1E-6,
        (sf_5e7.unloading_surface_pressure_mean + sf_5e7.unloading_surface_pressure_std) * 1E-6,
        color='C3', alpha=1, zorder=4)
    ax_particle_pressure.set_xlim(xmin=104.8, xmax=135)
    ax_particle_pressure.set_ylim(ymin=0, ymax=250)
    ax_particle_pressure.xaxis.set_major_locator(MultipleLocator(5))
    ax_particle_pressure.yaxis.set_major_locator(MultipleLocator(25))
    ax_particle_pressure.set_ylabel('Calendering surface pressure [MPa]')
    ax_particle_pressure.set_xlabel('Calendering surface height [µm]')
    fig_particle_pressure.tight_layout()
    ax_particle_pressure.legend(loc='best')
    fname = fig_dir + 'particle_pressure.svg'
    plt.savefig(fname)

    # =PARTICLE PRESSURE SPLIT==========================================================================================
    fig_particle_pressure_split, ax_particle_pressure_split = plt.subplots()
    ax_particle_pressure_split.set_xlabel('Calendering surface height [µm]')
    ax_particle_pressure_split.set_ylabel('Calendering surface pressure [MPa]')
    # ax_particle_pressure_split.set_xlim(xmin=0)
    # ax_particle_pressure_split.set_ylim(ymin=0, ymax=180)

    ax_particle_pressure_split.plot(ref.loading_surface_position, ref.loading_surface_pressure_matrix.T * 1E-6,
                                    color='C0', label=ref.label)
    ax_particle_pressure_split.plot(ref.unloading_surface_position,
                                    ref.unloading_surface_pressure_matrix.T * 1E-6, color='C0')

    ax_particle_pressure_split.plot(sf_5e7.loading_surface_position, sf_5e7.loading_surface_pressure_matrix.T * 1E-6,
                                    color='C1', label=sf_5e7.label)
    ax_particle_pressure_split.plot(sf_5e7.unloading_surface_position,
                                    sf_5e7.unloading_surface_pressure_matrix.T * 1E-6, color='C1')

    ax_particle_pressure_split.plot(sf_1e8.loading_surface_position, sf_1e8.loading_surface_pressure_matrix.T * 1E-6,
                                    color='C2', label=sf_1e8.label)
    ax_particle_pressure_split.plot(sf_1e8.unloading_surface_position,
                                    sf_1e8.unloading_surface_pressure_matrix.T * 1E-6, color='C2')

    ax_particle_pressure_split.plot(sf_2e8.loading_surface_position, sf_2e8.loading_surface_pressure_matrix.T * 1E-6,
                                    color='C3', label=sf_2e8.label)
    ax_particle_pressure_split.plot(sf_2e8.unloading_surface_position,
                                    sf_2e8.unloading_surface_pressure_matrix.T * 1E-6, color='C3')
    ax_particle_pressure_split.legend(loc='best')
    fig_particle_pressure_split.tight_layout()
    fname = fig_dir + 'particle_pressure_split.svg'
    plt.savefig(fname)

    # =PARTICLE CONTACTS_P ===============================================================================================
    fig_particle_contacts_p, ax_particle_contacts_p = plt.subplots()
    ax_particle_contacts_p.plot(ref.loading_surface_position * 1E2,
                                ref.loading_particle_contacts_mean,
                                'k-', zorder=1)
    ax_particle_contacts_p.plot(ref.unloading_surface_position * 1E2,
                                ref.unloading_particle_contacts_mean,
                                'k-', zorder=1)
    ax_particle_contacts_p.fill_between(
        ref.loading_surface_position * 1E2,
        (ref.loading_particle_contacts_mean - ref.loading_particle_contacts_std),
        (ref.loading_particle_contacts_mean + ref.loading_particle_contacts_std),
        color='C0', label=ref.label, alpha=1, zorder=1)
    ax_particle_contacts_p.fill_between(
        ref.unloading_surface_position * 1E2,
        (ref.unloading_particle_contacts_mean - ref.unloading_particle_contacts_std),
        (ref.unloading_particle_contacts_mean + ref.unloading_particle_contacts_std),
        color='C0', alpha=1, zorder=1)

    ax_particle_contacts_p.plot(sf_2e8.loading_surface_position * 1E2,
                                sf_2e8.loading_particle_contacts_mean,
                                'k-', zorder=2)
    ax_particle_contacts_p.plot(sf_2e8.unloading_surface_position * 1E2,
                                sf_2e8.unloading_particle_contacts_mean,
                                'k-', zorder=2)
    ax_particle_contacts_p.fill_between(
        sf_2e8.loading_surface_position * 1E2,
        (sf_2e8.loading_particle_contacts_mean - sf_2e8.loading_particle_contacts_std),
        (sf_2e8.loading_particle_contacts_mean + sf_2e8.loading_particle_contacts_std),
        color='C1', label=sf_2e8.label, alpha=1, zorder=2)
    ax_particle_contacts_p.fill_between(
        sf_2e8.unloading_surface_position * 1E2,
        (sf_2e8.unloading_particle_contacts_mean - sf_2e8.unloading_particle_contacts_std),
        (sf_2e8.unloading_particle_contacts_mean + sf_2e8.unloading_particle_contacts_std),
        color='C1', alpha=1, zorder=2)

    ax_particle_contacts_p.plot(sf_1e8.loading_surface_position * 1E2,
                                sf_1e8.loading_particle_contacts_mean,
                                'k-', zorder=3)
    ax_particle_contacts_p.plot(sf_1e8.unloading_surface_position * 1E2,
                                sf_1e8.unloading_particle_contacts_mean,
                                'k-', zorder=3)
    ax_particle_contacts_p.fill_between(
        sf_1e8.loading_surface_position * 1E2,
        (sf_1e8.loading_particle_contacts_mean - sf_1e8.loading_particle_contacts_std),
        (sf_1e8.loading_particle_contacts_mean + sf_1e8.loading_particle_contacts_std),
        color='C2', label=sf_1e8.label, alpha=1, zorder=3)
    ax_particle_contacts_p.fill_between(
        sf_1e8.unloading_surface_position * 1E2,
        (sf_1e8.unloading_particle_contacts_mean - sf_1e8.unloading_particle_contacts_std),
        (sf_1e8.unloading_particle_contacts_mean + sf_1e8.unloading_particle_contacts_std),
        color='C2', alpha=1, zorder=3)

    ax_particle_contacts_p.plot(sf_5e7.loading_surface_position * 1E2,
                                sf_5e7.loading_particle_contacts_mean,
                                'k-', zorder=4)
    ax_particle_contacts_p.plot(sf_5e7.unloading_surface_position * 1E2,
                                sf_5e7.unloading_particle_contacts_mean,
                                'k-', zorder=4)
    ax_particle_contacts_p.fill_between(
        sf_5e7.loading_surface_position * 1E2,
        (sf_5e7.loading_particle_contacts_mean - sf_5e7.loading_particle_contacts_std),
        (sf_5e7.loading_particle_contacts_mean + sf_5e7.loading_particle_contacts_std),
        color='C3', label=sf_5e7.label, alpha=1, zorder=4)
    ax_particle_contacts_p.fill_between(
        sf_5e7.unloading_surface_position * 1E2,
        (sf_5e7.unloading_particle_contacts_mean - sf_5e7.unloading_particle_contacts_std),
        (sf_5e7.unloading_particle_contacts_mean + sf_5e7.unloading_particle_contacts_std),
        color='C3', alpha=1, zorder=4)

    ax_particle_contacts_p.set_xlim(xmin=104.8, xmax=135)
    ax_particle_contacts_p.set_ylim(ymin=0, ymax=6000)
    ax_particle_contacts_p.xaxis.set_major_locator(MultipleLocator(5))
    ax_particle_contacts_p.yaxis.set_major_locator(MultipleLocator(500))
    ax_particle_contacts_p.set_ylabel('Particle contacts [-]')
    ax_particle_contacts_p.set_xlabel('Calendering surface height [µm]')
    ax_particle_contacts_p.legend(loc='best')
    fig_particle_contacts_p.tight_layout()
    fname = fig_dir + 'particle_contacts_p.svg'
    plt.savefig(fname)

    # =PARTICLE CONTACTS_P SPLIT========================================================================================
    fig_particle_contacts_p_split, ax_particle_contacts_p_split = plt.subplots()
    ax_particle_contacts_p_split.set_xlabel('Calendering surface height [µm]')
    ax_particle_contacts_p_split.set_ylabel('Particle contacts [-]')
    # ax_particle_contacts_p_split.set_xlim(xmin=0)
    # ax_particle_contacts_p_split.set_ylim(ymin=0, ymax=180)
    ax_particle_contacts_p_split.plot(ref.loading_surface_position, ref.loading_particle_contacts_matrix.T,
                                      color='C0', label=ref.label)
    ax_particle_contacts_p_split.plot(ref.unloading_surface_position,
                                      ref.unloading_particle_contacts_matrix.T, color='C0')

    ax_particle_contacts_p_split.plot(sf_5e7.loading_surface_position, sf_5e7.loading_particle_contacts_matrix.T,
                                      color='C1', label=sf_5e7.label)
    ax_particle_contacts_p_split.plot(sf_5e7.unloading_surface_position,
                                      sf_5e7.unloading_particle_contacts_matrix.T, color='C1')

    ax_particle_contacts_p_split.plot(sf_1e8.loading_surface_position, sf_1e8.loading_particle_contacts_matrix.T,
                                      color='C2', label=sf_1e8.label)
    ax_particle_contacts_p_split.plot(sf_1e8.unloading_surface_position,
                                      sf_1e8.unloading_particle_contacts_matrix.T, color='C2')

    ax_particle_contacts_p_split.plot(sf_2e8.loading_surface_position, sf_2e8.loading_particle_contacts_matrix.T,
                                      color='C3', label=sf_2e8.label)
    ax_particle_contacts_p_split.plot(sf_2e8.unloading_surface_position,
                                      sf_2e8.unloading_particle_contacts_matrix.T, color='C3')
    ax_particle_contacts_p_split.legend(loc='best')
    fig_particle_contacts_p_split.tight_layout()
    fname = fig_dir + 'particle_contacts_p_split.svg'
    plt.savefig(fname)

    # =PARTICLE CONTACTS_B==============================================================================================
    fig_particle_contacts_b, ax_particle_contacts_b = plt.subplots()
    ax_particle_contacts_b.plot(ref.loading_surface_position * 1E2,
                                ref.loading_binder_contacts_mean,
                                'k-', zorder=1)
    ax_particle_contacts_b.plot(ref.unloading_surface_position * 1E2,
                                ref.unloading_binder_contacts_mean,
                                'k-', zorder=1)
    ax_particle_contacts_b.fill_between(
        ref.loading_surface_position * 1E2,
        (ref.loading_binder_contacts_mean - ref.loading_binder_contacts_std),
        (ref.loading_binder_contacts_mean + ref.loading_binder_contacts_std),
        color='C0', label=ref.label, alpha=1, zorder=1)
    ax_particle_contacts_b.fill_between(
        ref.unloading_surface_position * 1E2,
        (ref.unloading_binder_contacts_mean - ref.unloading_binder_contacts_std),
        (ref.unloading_binder_contacts_mean + ref.unloading_binder_contacts_std),
        color='C0', alpha=1, zorder=1)

    ax_particle_contacts_b.plot(sf_2e8.loading_surface_position * 1E2,
                                sf_2e8.loading_binder_contacts_mean,
                                'k-', zorder=2)
    ax_particle_contacts_b.plot(sf_2e8.unloading_surface_position * 1E2,
                                sf_2e8.unloading_binder_contacts_mean,
                                'k-', zorder=2)
    ax_particle_contacts_b.fill_between(
        sf_2e8.loading_surface_position * 1E2,
        (sf_2e8.loading_binder_contacts_mean - sf_2e8.loading_binder_contacts_std),
        (sf_2e8.loading_binder_contacts_mean + sf_2e8.loading_binder_contacts_std),
        color='C1', label=sf_2e8.label, alpha=1, zorder=2)
    ax_particle_contacts_b.fill_between(
        sf_2e8.unloading_surface_position * 1E2,
        (sf_2e8.unloading_binder_contacts_mean - sf_2e8.unloading_binder_contacts_std),
        (sf_2e8.unloading_binder_contacts_mean + sf_2e8.unloading_binder_contacts_std),
        color='C1', alpha=1, zorder=2)

    ax_particle_contacts_b.plot(sf_1e8.loading_surface_position * 1E2,
                                sf_1e8.loading_binder_contacts_mean,
                                'k-', zorder=3)
    ax_particle_contacts_b.plot(sf_1e8.unloading_surface_position * 1E2,
                                sf_1e8.unloading_binder_contacts_mean,
                                'k-', zorder=3)
    ax_particle_contacts_b.fill_between(
        sf_1e8.loading_surface_position * 1E2,
        (sf_1e8.loading_binder_contacts_mean - sf_1e8.loading_binder_contacts_std),
        (sf_1e8.loading_binder_contacts_mean + sf_1e8.loading_binder_contacts_std),
        color='C2', label=sf_1e8.label, alpha=1, zorder=3)
    ax_particle_contacts_b.fill_between(
        sf_1e8.unloading_surface_position * 1E2,
        (sf_1e8.unloading_binder_contacts_mean - sf_1e8.unloading_binder_contacts_std),
        (sf_1e8.unloading_binder_contacts_mean + sf_1e8.unloading_binder_contacts_std),
        color='C2', alpha=1, zorder=3)

    ax_particle_contacts_b.plot(sf_5e7.loading_surface_position * 1E2,
                                sf_5e7.loading_binder_contacts_mean,
                                'k-', zorder=4)
    ax_particle_contacts_b.plot(sf_5e7.unloading_surface_position * 1E2,
                                sf_5e7.unloading_binder_contacts_mean,
                                'k-', zorder=4)
    ax_particle_contacts_b.fill_between(
        sf_5e7.loading_surface_position * 1E2,
        (sf_5e7.loading_binder_contacts_mean - sf_5e7.loading_binder_contacts_std),
        (sf_5e7.loading_binder_contacts_mean + sf_5e7.loading_binder_contacts_std),
        color='C3', label=sf_5e7.label, alpha=1, zorder=4)
    ax_particle_contacts_b.fill_between(
        sf_5e7.unloading_surface_position * 1E2,
        (sf_5e7.unloading_binder_contacts_mean - sf_5e7.unloading_binder_contacts_std),
        (sf_5e7.unloading_binder_contacts_mean + sf_5e7.unloading_binder_contacts_std),
        color='C3', alpha=1, zorder=4)

    ax_particle_contacts_b.set_xlim(xmin=104.8, xmax=135)
    ax_particle_contacts_b.set_ylim(ymin=8800, ymax=10200)
    ax_particle_contacts_b.xaxis.set_major_locator(MultipleLocator(5))
    ax_particle_contacts_b.yaxis.set_major_locator(MultipleLocator(200))
    ax_particle_contacts_b.set_ylabel('Binder contacts [-]')
    ax_particle_contacts_b.set_xlabel('Calendering surface height [µm]')
    fig_particle_contacts_b.tight_layout()
    ax_particle_contacts_b.legend(loc='best')
    fname = fig_dir + 'particle_contacts_b.svg'
    plt.savefig(fname)

    # =PARTICLE CONTACTS_B SPLIT========================================================================================
    fig_particle_contacts_b_split, ax_particle_contacts_b_split = plt.subplots()
    ax_particle_contacts_b_split.set_xlabel('Calendering surface height [µm]')
    ax_particle_contacts_b_split.set_ylabel('Binder contacts [-]')
    # ax_particle_contacts_b_split.set_xlim(xmin=0)
    # ax_particle_contacts_b_split.set_ylim(ymin=0, ymax=180)
    ax_particle_contacts_b_split.plot(ref.loading_surface_position, ref.loading_binder_contacts_matrix.T,
                                      color='C0', label=ref.label)
    ax_particle_contacts_b_split.plot(ref.unloading_surface_position,
                                      ref.unloading_binder_contacts_matrix.T, color='C0')

    ax_particle_contacts_b_split.plot(sf_5e7.loading_surface_position, sf_5e7.loading_binder_contacts_matrix.T,
                                      color='C1', label=sf_5e7.label)
    ax_particle_contacts_b_split.plot(sf_5e7.unloading_surface_position,
                                      sf_5e7.unloading_binder_contacts_matrix.T, color='C1')
    ax_particle_contacts_b_split.plot(sf_1e8.loading_surface_position, sf_1e8.loading_binder_contacts_matrix.T,
                                      color='C2', label=sf_1e8.label)
    ax_particle_contacts_b_split.plot(sf_1e8.unloading_surface_position,
                                      sf_1e8.unloading_binder_contacts_matrix.T, color='C2')

    ax_particle_contacts_b_split.plot(sf_2e8.loading_surface_position, sf_2e8.loading_binder_contacts_matrix.T,
                                      color='C3', label=sf_2e8.label)
    ax_particle_contacts_b_split.plot(sf_2e8.unloading_surface_position,
                                      sf_2e8.unloading_binder_contacts_matrix.T, color='C3')
    ax_particle_contacts_b_split.legend(loc='best')
    fig_particle_contacts_b_split.tight_layout()
    fname = fig_dir + 'particle_contacts_b_split.svg'
    plt.savefig(fname)

    # =PARTICLE FRACTURE_P ===============================================================================================
    fig_particle_fracture_p, ax_particle_fracture_p = plt.subplots()
    ax_particle_fracture_p.plot(sf_2e8.loading_surface_position * 1E2,
                                sf_2e8.loading_fractured_particles_mean * 1E2,
                                'k-', zorder=1)
    ax_particle_fracture_p.plot(sf_2e8.unloading_surface_position * 1E2,
                                sf_2e8.unloading_fractured_particles_mean * 1E2,
                                'k-', zorder=1)
    ax_particle_fracture_p.fill_between(
        sf_2e8.loading_surface_position * 1E2,
        (sf_2e8.loading_fractured_particles_mean - sf_2e8.loading_fractured_particles_std) * 1E2,
        (sf_2e8.loading_fractured_particles_mean + sf_2e8.loading_fractured_particles_std) * 1E2,
        color='C1', label=sf_2e8.label, alpha=1, zorder=1)
    ax_particle_fracture_p.fill_between(
        sf_2e8.unloading_surface_position * 1E2,
        (sf_2e8.unloading_fractured_particles_mean - sf_2e8.unloading_fractured_particles_std) * 1E2,
        (sf_2e8.unloading_fractured_particles_mean + sf_2e8.unloading_fractured_particles_std) * 1E2,
        color='C1', alpha=1, zorder=1)

    ax_particle_fracture_p.plot(sf_1e8.loading_surface_position * 1E2,
                                sf_1e8.loading_fractured_particles_mean * 1E2,
                                'k-', zorder=2)
    ax_particle_fracture_p.plot(sf_1e8.unloading_surface_position * 1E2,
                                sf_1e8.unloading_fractured_particles_mean * 1E2,
                                'k-', zorder=2)
    ax_particle_fracture_p.fill_between(
        sf_1e8.loading_surface_position * 1E2,
        (sf_1e8.loading_fractured_particles_mean - sf_1e8.loading_fractured_particles_std) * 1E2,
        (sf_1e8.loading_fractured_particles_mean + sf_1e8.loading_fractured_particles_std) * 1E2,
        color='C2', label=sf_1e8.label, alpha=1, zorder=2)
    ax_particle_fracture_p.fill_between(
        sf_1e8.unloading_surface_position * 1E2,
        (sf_1e8.unloading_fractured_particles_mean - sf_1e8.unloading_fractured_particles_std) * 1E2,
        (sf_1e8.unloading_fractured_particles_mean + sf_1e8.unloading_fractured_particles_std) * 1E2,
        color='C2', alpha=1, zorder=2)

    ax_particle_fracture_p.plot(sf_5e7.loading_surface_position * 1E2,
                                sf_5e7.loading_fractured_particles_mean * 1E2,
                                'k-', zorder=3)
    ax_particle_fracture_p.plot(sf_5e7.unloading_surface_position * 1E2,
                                sf_5e7.unloading_fractured_particles_mean * 1E2,
                                'k-', zorder=3)
    ax_particle_fracture_p.fill_between(
        sf_5e7.loading_surface_position * 1E2,
        (sf_5e7.loading_fractured_particles_mean - sf_5e7.loading_fractured_particles_std) * 1E2,
        (sf_5e7.loading_fractured_particles_mean + sf_5e7.loading_fractured_particles_std) * 1E2,
        color='C3', label=sf_5e7.label, alpha=1, zorder=3)
    ax_particle_fracture_p.fill_between(
        sf_5e7.unloading_surface_position * 1E2,
        (sf_5e7.unloading_fractured_particles_mean - sf_5e7.unloading_fractured_particles_std) * 1E2,
        (sf_5e7.unloading_fractured_particles_mean + sf_5e7.unloading_fractured_particles_std) * 1E2,
        color='C3', alpha=1, zorder=3)

    ax_particle_fracture_p.set_xlim(xmin=104.8, xmax=135)
    ax_particle_fracture_p.xaxis.set_major_locator(MultipleLocator(5))
    ax_particle_fracture_p.set_ylim(ymin=0, ymax=100)
    ax_particle_fracture_p.yaxis.set_major_locator(MultipleLocator(10))
    ax_particle_fracture_p.set_ylabel('Fractured particles [%]')
    ax_particle_fracture_p.set_xlabel('Calendering surface height [µm]')
    ax_particle_fracture_p.legend(loc='best')
    fig_particle_fracture_p.tight_layout()
    fname = fig_dir + 'particle_fracture_p.svg'
    plt.savefig(fname)

    # =PARTICLE FRACTURE_P SPLIT========================================================================================
    fig_particle_fracture_p_split, ax_particle_fracture_p_split = plt.subplots()
    ax_particle_fracture_p_split.set_xlabel('Calendering surface height [µm]')
    ax_particle_fracture_p_split.set_ylabel('Fraction of fractured particles [-]')
    # ax_particle_fracture_p_split.set_xlim(xmin=0)
    # ax_particle_contacts_p_split.set_ylim(ymin=0, ymax=180)
    ax_particle_fracture_p_split.plot(sf_5e7.loading_surface_position, sf_5e7.loading_fractured_particles_matrix.T,
                                      color='C1', label=sf_5e7.label)
    ax_particle_fracture_p_split.plot(sf_5e7.unloading_surface_position,
                                      sf_5e7.unloading_fractured_particles_matrix.T, color='C1')
    ax_particle_fracture_p_split.plot(sf_1e8.loading_surface_position, sf_1e8.loading_fractured_particles_matrix.T,
                                      color='C2', label=sf_1e8.label)
    ax_particle_fracture_p_split.plot(sf_1e8.unloading_surface_position,
                                      sf_1e8.unloading_fractured_particles_matrix.T, color='C2')

    ax_particle_fracture_p_split.plot(sf_2e8.loading_surface_position, sf_2e8.loading_fractured_particles_matrix.T,
                                      color='C3', label=sf_2e8.label)
    ax_particle_fracture_p_split.plot(sf_2e8.unloading_surface_position,
                                      sf_2e8.unloading_fractured_particles_matrix.T, color='C3')

    ax_particle_fracture_p_split.legend(loc='best')
    fig_particle_fracture_p_split.tight_layout()
    fname = fig_dir + 'particle_fracture_p_split.svg'
    plt.savefig(fname)

    # =PARTICLE FRACTURE_B==============================================================================================
    fig_particle_fracture_b, ax_particle_fracture_b = plt.subplots()
    ax_particle_fracture_b.plot(sf_2e8.loading_surface_position * 1E2,
                                sf_2e8.loading_fractured_binders_mean * 1E2,
                                'k-', zorder=1)
    ax_particle_fracture_b.plot(sf_2e8.unloading_surface_position * 1E2,
                                sf_2e8.unloading_fractured_binders_mean * 1E2,
                                'k-', zorder=1)
    ax_particle_fracture_b.fill_between(
        sf_2e8.loading_surface_position * 1E2,
        (sf_2e8.loading_fractured_binders_mean - sf_2e8.loading_fractured_binders_std) * 1E2,
        (sf_2e8.loading_fractured_binders_mean + sf_2e8.loading_fractured_binders_std) * 1E2,
        color='C1', label=sf_2e8.label, alpha=1, zorder=1)
    ax_particle_fracture_b.fill_between(
        sf_2e8.unloading_surface_position * 1E2,
        (sf_2e8.unloading_fractured_binders_mean - sf_2e8.unloading_fractured_binders_std) * 1E2,
        (sf_2e8.unloading_fractured_binders_mean + sf_2e8.unloading_fractured_binders_std) * 1E2,
        color='C1', alpha=1, zorder=1)

    ax_particle_fracture_b.plot(sf_1e8.loading_surface_position * 1E2,
                                sf_1e8.loading_fractured_binders_mean * 1E2,
                                'k-', zorder=2)
    ax_particle_fracture_b.plot(sf_1e8.unloading_surface_position * 1E2,
                                sf_1e8.unloading_fractured_binders_mean * 1E2,
                                'k-', zorder=2)
    ax_particle_fracture_b.fill_between(
        sf_1e8.loading_surface_position * 1E2,
        (sf_1e8.loading_fractured_binders_mean - sf_1e8.loading_fractured_binders_std) * 1E2,
        (sf_1e8.loading_fractured_binders_mean + sf_1e8.loading_fractured_binders_std) * 1E2,
        color='C2', label=sf_1e8.label, alpha=1, zorder=2)
    ax_particle_fracture_b.fill_between(
        sf_1e8.unloading_surface_position * 1E2,
        (sf_1e8.unloading_fractured_binders_mean - sf_1e8.unloading_fractured_binders_std) * 1E2,
        (sf_1e8.unloading_fractured_binders_mean + sf_1e8.unloading_fractured_binders_std) * 1E2,
        color='C2', alpha=1, zorder=2)

    ax_particle_fracture_b.plot(sf_5e7.loading_surface_position * 1E2,
                                sf_5e7.loading_fractured_binders_mean * 1E2,
                                'k-', zorder=3)
    ax_particle_fracture_b.plot(sf_5e7.unloading_surface_position * 1E2,
                                sf_5e7.unloading_fractured_binders_mean * 1E2,
                                'k-', zorder=3)
    ax_particle_fracture_b.fill_between(
        sf_5e7.loading_surface_position * 1E2,
        (sf_5e7.loading_fractured_binders_mean - sf_5e7.loading_fractured_binders_std) * 1E2,
        (sf_5e7.loading_fractured_binders_mean + sf_5e7.loading_fractured_binders_std) * 1E2,
        color='C3', label=sf_5e7.label, alpha=1, zorder=3)
    ax_particle_fracture_b.fill_between(
        sf_5e7.unloading_surface_position * 1E2,
        (sf_5e7.unloading_fractured_binders_mean - sf_5e7.unloading_fractured_binders_std) * 1E2,
        (sf_5e7.unloading_fractured_binders_mean + sf_5e7.unloading_fractured_binders_std) * 1E2,
        color='C3', alpha=1, zorder=3)

    ax_particle_fracture_b.set_xlim(xmin=104.8, xmax=135)
    ax_particle_fracture_b.xaxis.set_major_locator(MultipleLocator(5))
    ax_particle_fracture_b.set_ylim(ymin=0)
    ax_particle_fracture_b.set_ylabel('Fractured binder contacts [%]')
    ax_particle_fracture_b.set_xlabel('Calendering surface height [µm]')
    fig_particle_fracture_b.tight_layout()
    ax_particle_fracture_b.legend(loc='best')
    fname = fig_dir + 'particle_fracture_b.svg'
    plt.savefig(fname)

    # =PARTICLE FRACTURE_B SPLIT========================================================================================
    fig_particle_fracture_b_split, ax_particle_fracture_b_split = plt.subplots()
    ax_particle_fracture_b_split.set_xlabel('Calendering surface height [µm]')
    ax_particle_fracture_b_split.set_ylabel('Fraction of fractured binder contacts [-]')
    # ax_particle_fracture_b_split.set_xlim(xmin=0)
    # ax_particle_fracture_b_split.set_ylim(ymin=0, ymax=180)
    ax_particle_fracture_b_split.plot(sf_5e7.loading_surface_position, sf_5e7.loading_fractured_binders_matrix.T,
                                      color='C1', label=sf_5e7.label)
    ax_particle_fracture_b_split.plot(sf_5e7.unloading_surface_position,
                                      sf_5e7.unloading_fractured_binders_matrix.T, color='C1')

    ax_particle_fracture_b_split.plot(sf_1e8.loading_surface_position, sf_1e8.loading_fractured_binders_matrix.T,
                                      color='C2', label=sf_1e8.label)
    ax_particle_fracture_b_split.plot(sf_1e8.unloading_surface_position,
                                      sf_1e8.unloading_fractured_binders_matrix.T, color='C2')

    ax_particle_fracture_b_split.plot(sf_2e8.loading_surface_position, sf_2e8.loading_fractured_binders_matrix.T,
                                      color='C3', label=sf_2e8.label)
    ax_particle_fracture_b_split.plot(sf_2e8.unloading_surface_position,
                                      sf_2e8.unloading_fractured_binders_matrix.T, color='C3')

    fig_particle_fracture_b_split.tight_layout()
    fname = fig_dir + 'particle_fracture_b_split.svg'
    plt.savefig(fname)
    ########################################################################################################################

    ##########################################BINDERN FRACTURE##############################################################
    ef_5e_2 = Simulation('ef_5e-2', no_sims, r'$\varepsilon_{f}=5$ %')
    ef_1e_1 = Simulation('ef_1e-1', no_sims, r'$\varepsilon_{f}=10$ %')
    ef_2e_1 = Simulation('ef_2e-1', no_sims, r'$\varepsilon_{f}=20$ %')

    # =FIGURE BINDER PRESSURE===========================================================================================
    fig_binder_pressure, ax_binder_pressure = plt.subplots()
    ax_binder_pressure.plot(ref.loading_surface_position * 1E2,
                            ref.loading_surface_pressure_mean * 1E-6,
                            'k-', zorder=1)
    ax_binder_pressure.plot(ref.unloading_surface_position * 1E2,
                            ref.unloading_surface_pressure_mean * 1E-6,
                            'k-', zorder=1)
    ax_binder_pressure.fill_between(
        ref.loading_surface_position * 1E2,
        (ref.loading_surface_pressure_mean - ref.loading_surface_pressure_std) * 1E-6,
        (ref.loading_surface_pressure_mean + ref.loading_surface_pressure_std) * 1E-6,
        color='C0', label=ref.label, alpha=1, zorder=1)
    ax_binder_pressure.fill_between(
        ref.unloading_surface_position * 1E2,
        (ref.unloading_surface_pressure_mean - ref.unloading_surface_pressure_std) * 1E-6,
        (ref.unloading_surface_pressure_mean + ref.unloading_surface_pressure_std) * 1E-6,
        color='C0', alpha=1, zorder=1)

    ax_binder_pressure.plot(ef_2e_1.loading_surface_position * 1E2,
                            ef_2e_1.loading_surface_pressure_mean * 1E-6,
                            'k-', zorder=2)
    ax_binder_pressure.plot(ef_2e_1.unloading_surface_position * 1E2,
                            ef_2e_1.unloading_surface_pressure_mean * 1E-6,
                            'k-', zorder=2)
    ax_binder_pressure.fill_between(
        ef_2e_1.loading_surface_position * 1E2,
        (ef_2e_1.loading_surface_pressure_mean - ef_2e_1.loading_surface_pressure_std) * 1E-6,
        (ef_2e_1.loading_surface_pressure_mean + ef_2e_1.loading_surface_pressure_std) * 1E-6,
        color='C1', label=ef_2e_1.label, alpha=1, zorder=2)
    ax_binder_pressure.fill_between(
        ef_2e_1.unloading_surface_position * 1E2,
        (ef_2e_1.unloading_surface_pressure_mean - ef_2e_1.unloading_surface_pressure_std) * 1E-6,
        (ef_2e_1.unloading_surface_pressure_mean + ef_2e_1.unloading_surface_pressure_std) * 1E-6,
        color='C1', alpha=1, zorder=2)

    ax_binder_pressure.plot(ef_1e_1.loading_surface_position * 1E2,
                            ef_1e_1.loading_surface_pressure_mean * 1E-6,
                            'k-', zorder=3)
    ax_binder_pressure.plot(ef_1e_1.unloading_surface_position * 1E2,
                            ef_1e_1.unloading_surface_pressure_mean * 1E-6,
                            'k-', zorder=3)
    ax_binder_pressure.fill_between(
        ef_1e_1.loading_surface_position * 1E2,
        (ef_1e_1.loading_surface_pressure_mean - ef_1e_1.loading_surface_pressure_std) * 1E-6,
        (ef_1e_1.loading_surface_pressure_mean + ef_1e_1.loading_surface_pressure_std) * 1E-6,
        color='C2', label=ef_1e_1.label, alpha=1, zorder=3)
    ax_binder_pressure.fill_between(
        ef_1e_1.unloading_surface_position * 1E2,
        (ef_1e_1.unloading_surface_pressure_mean - ef_1e_1.unloading_surface_pressure_std) * 1E-6,
        (ef_1e_1.unloading_surface_pressure_mean + ef_1e_1.unloading_surface_pressure_std) * 1E-6,
        color='C2', alpha=1, zorder=3)

    ax_binder_pressure.plot(ef_5e_2.loading_surface_position * 1E2,
                            ef_5e_2.loading_surface_pressure_mean * 1E-6,
                            'k-', zorder=4)
    ax_binder_pressure.plot(ef_5e_2.unloading_surface_position * 1E2,
                            ef_5e_2.unloading_surface_pressure_mean * 1E-6,
                            'k-', zorder=4)
    ax_binder_pressure.fill_between(
        ef_5e_2.loading_surface_position * 1E2,
        (ef_5e_2.loading_surface_pressure_mean - ef_5e_2.loading_surface_pressure_std) * 1E-6,
        (ef_5e_2.loading_surface_pressure_mean + ef_5e_2.loading_surface_pressure_std) * 1E-6,
        color='C3', label=ef_5e_2.label, alpha=1, zorder=4)
    ax_binder_pressure.fill_between(
        ef_5e_2.unloading_surface_position * 1E2,
        (ef_5e_2.unloading_surface_pressure_mean - ef_5e_2.unloading_surface_pressure_std) * 1E-6,
        (ef_5e_2.unloading_surface_pressure_mean + ef_5e_2.unloading_surface_pressure_std) * 1E-6,
        color='C3', alpha=1, zorder=4)
    ax_binder_pressure.set_xlim(xmin=104.8, xmax=135)
    ax_binder_pressure.set_ylim(ymin=0, ymax=250)
    ax_binder_pressure.xaxis.set_major_locator(MultipleLocator(5))
    ax_binder_pressure.yaxis.set_major_locator(MultipleLocator(25))
    ax_binder_pressure.set_ylabel('Calendering surface pressure [MPa]')
    ax_binder_pressure.set_xlabel('Calendering surface height [µm]')
    fig_binder_pressure.tight_layout()
    ax_binder_pressure.legend(loc='best')
    fname = fig_dir + 'binder_pressure.svg'
    plt.savefig(fname)

    # =BINDER PRESSURE SPLIT============================================================================================
    fig_binder_pressure_split, ax_binder_pressure_split = plt.subplots()
    ax_binder_pressure_split.set_xlabel('Calendering surface height [µm]')
    ax_binder_pressure_split.set_ylabel('Calendering surface pressure [MPa]')
    # ax_binder_pressure_split.set_xlim(xmin=0)
    # ax_binder_pressure_split.set_ylim(ymin=0, ymax=180)

    ax_binder_pressure_split.plot(ref.loading_surface_position, ref.loading_surface_pressure_matrix.T * 1E-6,
                                  color='C0', label=ref.label)
    ax_binder_pressure_split.plot(ref.unloading_surface_position,
                                  ref.unloading_surface_pressure_matrix.T * 1E-6, color='C0')

    ax_binder_pressure_split.plot(ef_5e_2.loading_surface_position, ef_5e_2.loading_surface_pressure_matrix.T * 1E-6,
                                  color='C1', label=ef_5e_2.label)
    ax_binder_pressure_split.plot(ef_5e_2.unloading_surface_position,
                                  ef_5e_2.unloading_surface_pressure_matrix.T * 1E-6, color='C1')

    ax_binder_pressure_split.plot(ef_1e_1.loading_surface_position, ef_1e_1.loading_surface_pressure_matrix.T * 1E-6,
                                  color='C2', label=ef_1e_1.label)
    ax_binder_pressure_split.plot(ef_1e_1.unloading_surface_position,
                                  ef_1e_1.unloading_surface_pressure_matrix.T * 1E-6, color='C2')

    ax_binder_pressure_split.plot(ef_2e_1.loading_surface_position, ef_2e_1.loading_surface_pressure_matrix.T * 1E-6,
                                  color='C3', label=ef_2e_1.label)
    ax_binder_pressure_split.plot(ef_2e_1.unloading_surface_position,
                                  ef_2e_1.unloading_surface_pressure_matrix.T * 1E-6, color='C3')
    ax_binder_pressure_split.legend(loc='best')
    fig_binder_pressure_split.tight_layout()
    fname = fig_dir + 'binder_pressure_split.svg'
    plt.savefig(fname)

    # =BINDER CONTACTS_P ===============================================================================================
    fig_binder_contacts_p, ax_binder_contacts_p = plt.subplots()
    ax_binder_contacts_p.plot(ref.loading_surface_position * 1E2,
                              ref.loading_particle_contacts_mean,
                              'k-', zorder=1)
    ax_binder_contacts_p.plot(ref.unloading_surface_position * 1E2,
                              ref.unloading_particle_contacts_mean,
                              'k-', zorder=1)
    ax_binder_contacts_p.fill_between(
        ref.loading_surface_position * 1E2,
        (ref.loading_particle_contacts_mean - ref.loading_particle_contacts_std),
        (ref.loading_particle_contacts_mean + ref.loading_particle_contacts_std),
        color='C0', label=ref.label, alpha=1, zorder=1)
    ax_binder_contacts_p.fill_between(
        ref.unloading_surface_position * 1E2,
        (ref.unloading_particle_contacts_mean - ref.unloading_particle_contacts_std),
        (ref.unloading_particle_contacts_mean + ref.unloading_particle_contacts_std),
        color='C0', alpha=1, zorder=1)

    ax_binder_contacts_p.plot(ef_2e_1.loading_surface_position * 1E2,
                              ef_2e_1.loading_particle_contacts_mean,
                              'k-', zorder=2)
    ax_binder_contacts_p.plot(ef_2e_1.unloading_surface_position * 1E2,
                              ef_2e_1.unloading_particle_contacts_mean,
                              'k-', zorder=2)
    ax_binder_contacts_p.fill_between(
        ef_2e_1.loading_surface_position * 1E2,
        (ef_2e_1.loading_particle_contacts_mean - ef_2e_1.loading_particle_contacts_std),
        (ef_2e_1.loading_particle_contacts_mean + ef_2e_1.loading_particle_contacts_std),
        color='C1', label=ef_2e_1.label, alpha=1, zorder=2)
    ax_binder_contacts_p.fill_between(
        ef_2e_1.unloading_surface_position * 1E2,
        (ef_2e_1.unloading_particle_contacts_mean - ef_2e_1.unloading_particle_contacts_std),
        (ef_2e_1.unloading_particle_contacts_mean + ef_2e_1.unloading_particle_contacts_std),
        color='C1', alpha=1, zorder=2)

    ax_binder_contacts_p.plot(ef_1e_1.loading_surface_position * 1E2,
                              ef_1e_1.loading_particle_contacts_mean,
                              'k-', zorder=3)
    ax_binder_contacts_p.plot(ef_1e_1.unloading_surface_position * 1E2,
                              ef_1e_1.unloading_particle_contacts_mean,
                              'k-', zorder=3)
    ax_binder_contacts_p.fill_between(
        ef_1e_1.loading_surface_position * 1E2,
        (ef_1e_1.loading_particle_contacts_mean - ef_1e_1.loading_particle_contacts_std),
        (ef_1e_1.loading_particle_contacts_mean + ef_1e_1.loading_particle_contacts_std),
        color='C2', label=ef_1e_1.label, alpha=1, zorder=3)
    ax_binder_contacts_p.fill_between(
        ef_1e_1.unloading_surface_position * 1E2,
        (ef_1e_1.unloading_particle_contacts_mean - ef_1e_1.unloading_particle_contacts_std),
        (ef_1e_1.unloading_particle_contacts_mean + ef_1e_1.unloading_particle_contacts_std),
        color='C2', alpha=1, zorder=3)

    ax_binder_contacts_p.plot(ef_5e_2.loading_surface_position * 1E2,
                              ef_5e_2.loading_particle_contacts_mean,
                              'k-', zorder=4)
    ax_binder_contacts_p.plot(ef_5e_2.unloading_surface_position * 1E2,
                              ef_5e_2.unloading_particle_contacts_mean,
                              'k-', zorder=4)
    ax_binder_contacts_p.fill_between(
        ef_5e_2.loading_surface_position * 1E2,
        (ef_5e_2.loading_particle_contacts_mean - ef_5e_2.loading_particle_contacts_std),
        (ef_5e_2.loading_particle_contacts_mean + ef_5e_2.loading_particle_contacts_std),
        color='C3', label=ef_5e_2.label, alpha=1, zorder=4)
    ax_binder_contacts_p.fill_between(
        ef_5e_2.unloading_surface_position * 1E2,
        (ef_5e_2.unloading_particle_contacts_mean - ef_5e_2.unloading_particle_contacts_std),
        (ef_5e_2.unloading_particle_contacts_mean + ef_5e_2.unloading_particle_contacts_std),
        color='C3', alpha=1, zorder=4)

    ax_binder_contacts_p.set_xlim(xmin=104.8, xmax=135)
    ax_binder_contacts_p.set_ylim(ymin=0, ymax=6000)
    ax_binder_contacts_p.xaxis.set_major_locator(MultipleLocator(5))
    ax_binder_contacts_p.yaxis.set_major_locator(MultipleLocator(500))
    ax_binder_contacts_p.set_ylabel('Particle contacts [-]')
    ax_binder_contacts_p.set_xlabel('Calendering surface height [µm]')
    ax_binder_contacts_p.legend(loc='best')
    fig_binder_contacts_p.tight_layout()
    fname = fig_dir + 'binder_contacts_p.svg'
    plt.savefig(fname)

    # =BINDER CONTACTS_P SPLIT==========================================================================================
    fig_binder_contacts_p_split, ax_binder_contacts_p_split = plt.subplots()
    ax_binder_contacts_p_split.set_xlabel('Calendering surface height [µm]')
    ax_binder_contacts_p_split.set_ylabel('Particle contacts [-]')
    # ax_binder_contacts_p_split.set_xlim(xmin=0)
    # ax_binder_contacts_p_split.set_ylim(ymin=0, ymax=180)
    ax_binder_contacts_p_split.plot(ref.loading_surface_position, ref.loading_particle_contacts_matrix.T,
                                    color='C0', label=ref.label)
    ax_binder_contacts_p_split.plot(ref.unloading_surface_position,
                                    ref.unloading_particle_contacts_matrix.T, color='C0')

    ax_binder_contacts_p_split.plot(ef_5e_2.loading_surface_position, ef_5e_2.loading_particle_contacts_matrix.T,
                                    color='C1', label=ef_5e_2.label)
    ax_binder_contacts_p_split.plot(ef_5e_2.unloading_surface_position,
                                    ef_5e_2.unloading_particle_contacts_matrix.T, color='C1')

    ax_binder_contacts_p_split.plot(ef_1e_1.loading_surface_position, ef_1e_1.loading_particle_contacts_matrix.T,
                                    color='C2', label=ef_1e_1.label)
    ax_binder_contacts_p_split.plot(ef_1e_1.unloading_surface_position,
                                    ef_1e_1.unloading_particle_contacts_matrix.T, color='C2')

    ax_binder_contacts_p_split.plot(ef_2e_1.loading_surface_position, ef_2e_1.loading_particle_contacts_matrix.T,
                                    color='C3', label=ef_2e_1.label)
    ax_binder_contacts_p_split.plot(ef_2e_1.unloading_surface_position,
                                    ef_2e_1.unloading_particle_contacts_matrix.T, color='C3')
    ax_binder_contacts_p_split.legend(loc='best')
    fig_binder_contacts_p_split.tight_layout()
    fname = fig_dir + 'binder_contacts_p_split.svg'
    plt.savefig(fname)

    # =BINDER CONTACTS_B================================================================================================
    fig_binder_contacts_b, ax_binder_contacts_b = plt.subplots()
    ax_binder_contacts_b.plot(ref.loading_surface_position * 1E2,
                              ref.loading_binder_contacts_mean,
                              'k-', zorder=1)
    ax_binder_contacts_b.plot(ref.unloading_surface_position * 1E2,
                              ref.unloading_binder_contacts_mean,
                              'k-', zorder=1)
    ax_binder_contacts_b.fill_between(
        ref.loading_surface_position * 1E2,
        (ref.loading_binder_contacts_mean - ref.loading_binder_contacts_std),
        (ref.loading_binder_contacts_mean + ref.loading_binder_contacts_std),
        color='C0', label=ref.label, alpha=1, zorder=1)
    ax_binder_contacts_b.fill_between(
        ref.unloading_surface_position * 1E2,
        (ref.unloading_binder_contacts_mean - ref.unloading_binder_contacts_std),
        (ref.unloading_binder_contacts_mean + ref.unloading_binder_contacts_std),
        color='C0', alpha=1, zorder=1)

    ax_binder_contacts_b.plot(ef_2e_1.loading_surface_position * 1E2,
                              ef_2e_1.loading_binder_contacts_mean,
                              'k-', zorder=2)
    ax_binder_contacts_b.plot(ef_2e_1.unloading_surface_position * 1E2,
                              ef_2e_1.unloading_binder_contacts_mean,
                              'k-', zorder=2)
    ax_binder_contacts_b.fill_between(
        ef_2e_1.loading_surface_position * 1E2,
        (ef_2e_1.loading_binder_contacts_mean - ef_2e_1.loading_binder_contacts_std),
        (ef_2e_1.loading_binder_contacts_mean + ef_2e_1.loading_binder_contacts_std),
        color='C1', label=ef_2e_1.label, alpha=1, zorder=2)
    ax_binder_contacts_b.fill_between(
        ef_2e_1.unloading_surface_position * 1E2,
        (ef_2e_1.unloading_binder_contacts_mean - ef_2e_1.unloading_binder_contacts_std),
        (ef_2e_1.unloading_binder_contacts_mean + ef_2e_1.unloading_binder_contacts_std),
        color='C1', alpha=1, zorder=2)

    ax_binder_contacts_b.plot(ef_1e_1.loading_surface_position * 1E2,
                              ef_1e_1.loading_binder_contacts_mean,
                              'k-', zorder=3)
    ax_binder_contacts_b.plot(ef_1e_1.unloading_surface_position * 1E2,
                              ef_1e_1.unloading_binder_contacts_mean,
                              'k-', zorder=3)
    ax_binder_contacts_b.fill_between(
        ef_1e_1.loading_surface_position * 1E2,
        (ef_1e_1.loading_binder_contacts_mean - ef_1e_1.loading_binder_contacts_std),
        (ef_1e_1.loading_binder_contacts_mean + ef_1e_1.loading_binder_contacts_std),
        color='C2', label=ef_1e_1.label, alpha=1, zorder=3)
    ax_binder_contacts_b.fill_between(
        ef_1e_1.unloading_surface_position * 1E2,
        (ef_1e_1.unloading_binder_contacts_mean - ef_1e_1.unloading_binder_contacts_std),
        (ef_1e_1.unloading_binder_contacts_mean + ef_1e_1.unloading_binder_contacts_std),
        color='C2', alpha=1, zorder=3)

    ax_binder_contacts_b.plot(ef_5e_2.loading_surface_position * 1E2,
                              ef_5e_2.loading_binder_contacts_mean,
                              'k-', zorder=4)
    ax_binder_contacts_b.plot(ef_5e_2.unloading_surface_position * 1E2,
                              ef_5e_2.unloading_binder_contacts_mean,
                              'k-', zorder=4)
    ax_binder_contacts_b.fill_between(
        ef_5e_2.loading_surface_position * 1E2,
        (ef_5e_2.loading_binder_contacts_mean - ef_5e_2.loading_binder_contacts_std),
        (ef_5e_2.loading_binder_contacts_mean + ef_5e_2.loading_binder_contacts_std),
        color='C3', label=ef_5e_2.label, alpha=1, zorder=4)
    ax_binder_contacts_b.fill_between(
        ef_5e_2.unloading_surface_position * 1E2,
        (ef_5e_2.unloading_binder_contacts_mean - ef_5e_2.unloading_binder_contacts_std),
        (ef_5e_2.unloading_binder_contacts_mean + ef_5e_2.unloading_binder_contacts_std),
        color='C3', alpha=1, zorder=4)

    ax_binder_contacts_b.set_xlim(xmin=104.8, xmax=135)
    ax_binder_contacts_b.set_ylim(ymin=8800, ymax=10200)
    ax_binder_contacts_b.xaxis.set_major_locator(MultipleLocator(5))
    ax_binder_contacts_b.yaxis.set_major_locator(MultipleLocator(200))
    ax_binder_contacts_b.set_ylabel('Binder contacts [-]')
    ax_binder_contacts_b.set_xlabel('Calendering surface height [µm]')
    fig_binder_contacts_b.tight_layout()
    ax_binder_contacts_b.legend(loc='best')
    fname = fig_dir + 'binder_contacts_b.svg'
    plt.savefig(fname)

    # =BINDER CONTACTS_B SPLIT==========================================================================================
    fig_binder_contacts_b_split, ax_binder_contacts_b_split = plt.subplots()
    ax_binder_contacts_b_split.set_xlabel('Calendering surface height [µm]')
    ax_binder_contacts_b_split.set_ylabel('Binder contacts [-]')
    # ax_binder_contacts_b_split.set_xlim(xmin=0)
    # ax_binder_contacts_b_split.set_ylim(ymin=0, ymax=180)
    ax_binder_contacts_b_split.plot(ref.loading_surface_position, ref.loading_binder_contacts_matrix.T,
                                    color='C0', label=ref.label)
    ax_binder_contacts_b_split.plot(ref.unloading_surface_position,
                                    ref.unloading_binder_contacts_matrix.T, color='C0')

    ax_binder_contacts_b_split.plot(ef_5e_2.loading_surface_position, ef_5e_2.loading_binder_contacts_matrix.T,
                                    color='C1', label=ef_5e_2.label)
    ax_binder_contacts_b_split.plot(ef_5e_2.unloading_surface_position,
                                    ef_5e_2.unloading_binder_contacts_matrix.T, color='C1')
    ax_binder_contacts_b_split.plot(ef_1e_1.loading_surface_position, ef_1e_1.loading_binder_contacts_matrix.T,
                                    color='C2', label=ef_1e_1.label)
    ax_binder_contacts_b_split.plot(ef_1e_1.unloading_surface_position,
                                    ef_1e_1.unloading_binder_contacts_matrix.T, color='C2')

    ax_binder_contacts_b_split.plot(ef_2e_1.loading_surface_position, ef_2e_1.loading_binder_contacts_matrix.T,
                                    color='C3', label=ef_2e_1.label)
    ax_binder_contacts_b_split.plot(ef_2e_1.unloading_surface_position,
                                    ef_2e_1.unloading_binder_contacts_matrix.T, color='C3')
    ax_binder_contacts_b_split.legend(loc='best')
    fig_binder_contacts_b_split.tight_layout()
    fname = fig_dir + 'binder_contacts_b_split.svg'
    plt.savefig(fname)

    # =PARTICLE FRACTURE_P ===============================================================================================
    fig_binder_fracture_p, ax_binder_fracture_p = plt.subplots()
    ax_binder_fracture_p.plot(ef_2e_1.loading_surface_position * 1E2,
                              ef_2e_1.loading_fractured_particles_mean * 1E2,
                              'k-', zorder=1)
    ax_binder_fracture_p.plot(ef_2e_1.unloading_surface_position * 1E2,
                              ef_2e_1.unloading_fractured_particles_mean * 1E2,
                              'k-', zorder=1)
    ax_binder_fracture_p.fill_between(
        ef_2e_1.loading_surface_position * 1E2,
        (ef_2e_1.loading_fractured_particles_mean - ef_2e_1.loading_fractured_particles_std) * 1E2,
        (ef_2e_1.loading_fractured_particles_mean + ef_2e_1.loading_fractured_particles_std) * 1E2,
        color='C1', label=ef_2e_1.label, alpha=1, zorder=1)
    ax_binder_fracture_p.fill_between(
        ef_2e_1.unloading_surface_position * 1E2,
        (ef_2e_1.unloading_fractured_particles_mean - ef_2e_1.unloading_fractured_particles_std) * 1E2,
        (ef_2e_1.unloading_fractured_particles_mean + ef_2e_1.unloading_fractured_particles_std) * 1E2,
        color='C1', alpha=1, zorder=1)

    ax_binder_fracture_p.plot(ef_1e_1.loading_surface_position * 1E2,
                              ef_1e_1.loading_fractured_particles_mean * 1E2,
                              'k-', zorder=2)
    ax_binder_fracture_p.plot(ef_1e_1.unloading_surface_position * 1E2,
                              ef_1e_1.unloading_fractured_particles_mean * 1E2,
                              'k-', zorder=2)
    ax_binder_fracture_p.fill_between(
        ef_1e_1.loading_surface_position * 1E2,
        (ef_1e_1.loading_fractured_particles_mean - ef_1e_1.loading_fractured_particles_std) * 1E2,
        (ef_1e_1.loading_fractured_particles_mean + ef_1e_1.loading_fractured_particles_std) * 1E2,
        color='C2', label=ef_1e_1.label, alpha=1, zorder=2)
    ax_binder_fracture_p.fill_between(
        ef_1e_1.unloading_surface_position * 1E2,
        (ef_1e_1.unloading_fractured_particles_mean - ef_1e_1.unloading_fractured_particles_std) * 1E2,
        (ef_1e_1.unloading_fractured_particles_mean + ef_1e_1.unloading_fractured_particles_std) * 1E2,
        color='C2', alpha=1, zorder=2)

    ax_binder_fracture_p.plot(ef_5e_2.loading_surface_position * 1E2,
                              ef_5e_2.loading_fractured_particles_mean * 1E2,
                              'k-', zorder=3)
    ax_binder_fracture_p.plot(ef_5e_2.unloading_surface_position * 1E2,
                              ef_5e_2.unloading_fractured_particles_mean * 1E2,
                              'k-', zorder=3)
    ax_binder_fracture_p.fill_between(
        ef_5e_2.loading_surface_position * 1E2,
        (ef_5e_2.loading_fractured_particles_mean - ef_5e_2.loading_fractured_particles_std) * 1E2,
        (ef_5e_2.loading_fractured_particles_mean + ef_5e_2.loading_fractured_particles_std) * 1E2,
        color='C3', label=ef_5e_2.label, alpha=1, zorder=3)
    ax_binder_fracture_p.fill_between(
        ef_5e_2.unloading_surface_position * 1E2,
        (ef_5e_2.unloading_fractured_particles_mean - ef_5e_2.unloading_fractured_particles_std) * 1E2,
        (ef_5e_2.unloading_fractured_particles_mean + ef_5e_2.unloading_fractured_particles_std) * 1E2,
        color='C3', alpha=1, zorder=3)

    ax_binder_fracture_p.set_xlim(xmin=104.8, xmax=135)
    ax_binder_fracture_p.xaxis.set_major_locator(MultipleLocator(5))
    ax_binder_fracture_p.set_ylim(ymin=0, ymax=100)
    ax_binder_fracture_p.yaxis.set_major_locator(MultipleLocator(10))
    ax_binder_fracture_p.set_ylabel('Fractured particles [%]')
    ax_binder_fracture_p.set_xlabel('Calendering surface height [µm]')
    ax_binder_fracture_p.legend(loc='best')
    fig_binder_fracture_p.tight_layout()
    fname = fig_dir + 'binder_fracture_p.svg'
    plt.savefig(fname)

    # =BINDER FRACTURE_P SPLIT========================================================================================
    fig_binder_fracture_p_split, ax_binder_fracture_p_split = plt.subplots()
    ax_binder_fracture_p_split.set_xlabel('Calendering surface height [µm]')
    ax_binder_fracture_p_split.set_ylabel('Fraction of fractured particles [-]')
    # ax_binder_fracture_p_split.set_xlim(xmin=0)
    # ax_binder_contacts_p_split.set_ylim(ymin=0, ymax=180)
    ax_binder_fracture_p_split.plot(ef_5e_2.loading_surface_position, ef_5e_2.loading_fractured_particles_matrix.T,
                                    color='C1', label=ef_5e_2.label)
    ax_binder_fracture_p_split.plot(ef_5e_2.unloading_surface_position,
                                    ef_5e_2.unloading_fractured_particles_matrix.T, color='C1')
    ax_binder_fracture_p_split.plot(ef_1e_1.loading_surface_position, ef_1e_1.loading_fractured_particles_matrix.T,
                                    color='C2', label=ef_1e_1.label)
    ax_binder_fracture_p_split.plot(ef_1e_1.unloading_surface_position,
                                    ef_1e_1.unloading_fractured_particles_matrix.T, color='C2')

    ax_binder_fracture_p_split.plot(ef_2e_1.loading_surface_position, ef_2e_1.loading_fractured_particles_matrix.T,
                                    color='C3', label=ef_2e_1.label)
    ax_binder_fracture_p_split.plot(ef_2e_1.unloading_surface_position,
                                    ef_2e_1.unloading_fractured_particles_matrix.T, color='C3')

    ax_binder_fracture_p_split.legend(loc='best')
    fig_binder_fracture_p_split.tight_layout()
    fname = fig_dir + 'binder_fracture_p_split.svg'
    plt.savefig(fname)

    # =PARTICLE FRACTURE_B==============================================================================================
    fig_binder_fracture_b, ax_binder_fracture_b = plt.subplots()
    ax_binder_fracture_b.plot(ef_2e_1.loading_surface_position * 1E2,
                              ef_2e_1.loading_fractured_binders_mean * 1E2,
                              'k-', zorder=1)
    ax_binder_fracture_b.plot(ef_2e_1.unloading_surface_position * 1E2,
                              ef_2e_1.unloading_fractured_binders_mean * 1E2,
                              'k-', zorder=1)
    ax_binder_fracture_b.fill_between(
        ef_2e_1.loading_surface_position * 1E2,
        (ef_2e_1.loading_fractured_binders_mean - ef_2e_1.loading_fractured_binders_std) * 1E2,
        (ef_2e_1.loading_fractured_binders_mean + ef_2e_1.loading_fractured_binders_std) * 1E2,
        color='C1', label=ef_2e_1.label, alpha=1, zorder=1)
    ax_binder_fracture_b.fill_between(
        ef_2e_1.unloading_surface_position * 1E2,
        (ef_2e_1.unloading_fractured_binders_mean - ef_2e_1.unloading_fractured_binders_std) * 1E2,
        (ef_2e_1.unloading_fractured_binders_mean + ef_2e_1.unloading_fractured_binders_std) * 1E2,
        color='C1', alpha=1, zorder=1)

    ax_binder_fracture_b.plot(ef_1e_1.loading_surface_position * 1E2,
                              ef_1e_1.loading_fractured_binders_mean * 1E2,
                              'k-', zorder=2)
    ax_binder_fracture_b.plot(ef_1e_1.unloading_surface_position * 1E2,
                              ef_1e_1.unloading_fractured_binders_mean * 1E2,
                              'k-', zorder=2)
    ax_binder_fracture_b.fill_between(
        ef_1e_1.loading_surface_position * 1E2,
        (ef_1e_1.loading_fractured_binders_mean - ef_1e_1.loading_fractured_binders_std) * 1E2,
        (ef_1e_1.loading_fractured_binders_mean + ef_1e_1.loading_fractured_binders_std) * 1E2,
        color='C2', label=ef_1e_1.label, alpha=1, zorder=2)
    ax_binder_fracture_b.fill_between(
        ef_1e_1.unloading_surface_position * 1E2,
        (ef_1e_1.unloading_fractured_binders_mean - ef_1e_1.unloading_fractured_binders_std) * 1E2,
        (ef_1e_1.unloading_fractured_binders_mean + ef_1e_1.unloading_fractured_binders_std) * 1E2,
        color='C2', alpha=1, zorder=2)

    ax_binder_fracture_b.plot(ef_5e_2.loading_surface_position * 1E2,
                              ef_5e_2.loading_fractured_binders_mean * 1E2,
                              'k-', zorder=3)
    ax_binder_fracture_b.plot(ef_5e_2.unloading_surface_position * 1E2,
                              ef_5e_2.unloading_fractured_binders_mean * 1E2,
                              'k-', zorder=3)
    ax_binder_fracture_b.fill_between(
        ef_5e_2.loading_surface_position * 1E2,
        (ef_5e_2.loading_fractured_binders_mean - ef_5e_2.loading_fractured_binders_std) * 1E2,
        (ef_5e_2.loading_fractured_binders_mean + ef_5e_2.loading_fractured_binders_std) * 1E2,
        color='C3', label=ef_5e_2.label, alpha=1, zorder=3)
    ax_binder_fracture_b.fill_between(
        ef_5e_2.unloading_surface_position * 1E2,
        (ef_5e_2.unloading_fractured_binders_mean - ef_5e_2.unloading_fractured_binders_std) * 1E2,
        (ef_5e_2.unloading_fractured_binders_mean + ef_5e_2.unloading_fractured_binders_std) * 1E2,
        color='C3', alpha=1, zorder=3)

    ax_binder_fracture_b.set_xlim(xmin=104.8, xmax=135)
    ax_binder_fracture_b.xaxis.set_major_locator(MultipleLocator(5))
    ax_binder_fracture_b.set_ylim(ymin=0, ymax=25)
    ax_binder_fracture_b.yaxis.set_major_locator(MultipleLocator(5))
    ax_binder_fracture_b.set_ylabel('Fractured binder contacts [%]')
    ax_binder_fracture_b.set_xlabel('Calendering surface height [µm]')
    fig_binder_fracture_b.tight_layout()
    ax_binder_fracture_b.legend(loc='best')
    fname = fig_dir + 'binder_fracture_b.svg'
    plt.savefig(fname)

    # =PARTICLE FRACTURE_B SPLIT========================================================================================
    fig_binder_fracture_b_split, ax_binder_fracture_b_split = plt.subplots()
    ax_binder_fracture_b_split.set_xlabel('Calendering surface height [µm]')
    ax_binder_fracture_b_split.set_ylabel('Fraction of fractured binder contacts [-]')
    # ax_binder_fracture_b_split.set_xlim(xmin=0)
    # ax_binder_fracture_b_split.set_ylim(ymin=0, ymax=180)
    ax_binder_fracture_b_split.plot(ef_5e_2.loading_surface_position, ef_5e_2.loading_fractured_binders_matrix.T,
                                    color='C1', label=ef_5e_2.label)
    ax_binder_fracture_b_split.plot(ef_5e_2.unloading_surface_position,
                                    ef_5e_2.unloading_fractured_binders_matrix.T, color='C1')

    ax_binder_fracture_b_split.plot(ef_1e_1.loading_surface_position, ef_1e_1.loading_fractured_binders_matrix.T,
                                    color='C2', label=ef_1e_1.label)
    ax_binder_fracture_b_split.plot(ef_1e_1.unloading_surface_position,
                                    ef_1e_1.unloading_fractured_binders_matrix.T, color='C2')

    ax_binder_fracture_b_split.plot(ef_2e_1.loading_surface_position, ef_2e_1.loading_fractured_binders_matrix.T,
                                    color='C3', label=ef_2e_1.label)
    ax_binder_fracture_b_split.plot(ef_2e_1.unloading_surface_position,
                                    ef_2e_1.unloading_fractured_binders_matrix.T, color='C3')

    fig_binder_fracture_b_split.tight_layout()
    fname = fig_dir + 'binder_fracture_b_split.svg'
    plt.savefig(fname)
    ########################################################################################################################

    #######################################BINDER AND PARTICLE FRACTURE#####################################################

    sf_2e8_ef_1e_1 = Simulation('sf_2e8_ef_1e-1', no_sims, r'$\sigma_{f}=200$ MPa, $\varepsilon_{f}=10$ %')

    # =FIGURE BINDER PARTICLE PRESSURE==================================================================================
    fig_particle_binder_pressure, ax_particle_binder_pressure = plt.subplots()
    ax_particle_binder_pressure.plot(ref.loading_surface_position * 1E2,
                                     ref.loading_surface_pressure_mean * 1E-6,
                                     'k-', zorder=1)
    ax_particle_binder_pressure.plot(ref.unloading_surface_position * 1E2,
                                     ref.unloading_surface_pressure_mean * 1E-6,
                                     'k-', zorder=1)
    ax_particle_binder_pressure.fill_between(
        ref.loading_surface_position * 1E2,
        (ref.loading_surface_pressure_mean - ref.loading_surface_pressure_std) * 1E-6,
        (ref.loading_surface_pressure_mean + ref.loading_surface_pressure_std) * 1E-6,
        color='C0', label=ref.label, alpha=1, zorder=1)
    ax_particle_binder_pressure.fill_between(
        ref.unloading_surface_position * 1E2,
        (ref.unloading_surface_pressure_mean - ref.unloading_surface_pressure_std) * 1E-6,
        (ref.unloading_surface_pressure_mean + ref.unloading_surface_pressure_std) * 1E-6,
        color='C0', alpha=1, zorder=1)

    ax_particle_binder_pressure.plot(sf_2e8_ef_1e_1.loading_surface_position * 1E2,
                                     sf_2e8_ef_1e_1.loading_surface_pressure_mean * 1E-6,
                                     'k-', zorder=2)
    ax_particle_binder_pressure.plot(sf_2e8_ef_1e_1.unloading_surface_position * 1E2,
                                     sf_2e8_ef_1e_1.unloading_surface_pressure_mean * 1E-6,
                                     'k-', zorder=2)
    ax_particle_binder_pressure.fill_between(
        sf_2e8_ef_1e_1.loading_surface_position * 1E2,
        (sf_2e8_ef_1e_1.loading_surface_pressure_mean - sf_2e8_ef_1e_1.loading_surface_pressure_std) * 1E-6,
        (sf_2e8_ef_1e_1.loading_surface_pressure_mean + sf_2e8_ef_1e_1.loading_surface_pressure_std) * 1E-6,
        color='C3', label=sf_2e8_ef_1e_1.label, alpha=1, zorder=2)
    ax_particle_binder_pressure.fill_between(
        sf_2e8_ef_1e_1.unloading_surface_position * 1E2,
        (sf_2e8_ef_1e_1.unloading_surface_pressure_mean - sf_2e8_ef_1e_1.unloading_surface_pressure_std) * 1E-6,
        (sf_2e8_ef_1e_1.unloading_surface_pressure_mean + sf_2e8_ef_1e_1.unloading_surface_pressure_std) * 1E-6,
        color='C3', alpha=1, zorder=2)

    ax_particle_binder_pressure.set_xlim(xmin=104.8, xmax=135)
    ax_particle_binder_pressure.set_ylim(ymin=0, ymax=250)
    ax_particle_binder_pressure.xaxis.set_major_locator(MultipleLocator(5))
    ax_particle_binder_pressure.yaxis.set_major_locator(MultipleLocator(25))
    ax_particle_binder_pressure.set_ylabel('Calendering surface pressure [MPa]')
    ax_particle_binder_pressure.set_xlabel('Calendering surface height [µm]')
    fig_particle_binder_pressure.tight_layout()
    ax_particle_binder_pressure.legend(loc='best')
    fname = fig_dir + 'particle_binder_pressure.svg'
    plt.savefig(fname)

    # =BINDER PARTICLE PRESSURE SPLIT===================================================================================
    fig_particle_binder_pressure_split, ax_particle_binder_pressure_split = plt.subplots()
    ax_particle_binder_pressure_split.set_xlabel('Calendering surface height [µm]')
    ax_particle_binder_pressure_split.set_ylabel('Calendering surface pressure [MPa]')
    # ax_particle_binder_pressure_split.set_xlim(xmin=0)
    # ax_particle_binder_pressure_split.set_ylim(ymin=0, ymax=180)

    ax_particle_binder_pressure_split.plot(ref.loading_surface_position, ref.loading_surface_pressure_matrix.T * 1E-6,
                                           color='C0', label=ref.label)
    ax_particle_binder_pressure_split.plot(ref.unloading_surface_position,
                                           ref.unloading_surface_pressure_matrix.T * 1E-6, color='C0')

    ax_particle_binder_pressure_split.plot(sf_2e8_ef_1e_1.loading_surface_position,
                                           sf_2e8_ef_1e_1.loading_surface_pressure_matrix.T * 1E-6,
                                           color='C3', label=sf_2e8_ef_1e_1.label)
    ax_particle_binder_pressure_split.plot(sf_2e8_ef_1e_1.unloading_surface_position,
                                           sf_2e8_ef_1e_1.unloading_surface_pressure_matrix.T * 1E-6, color='C3')

    ax_particle_binder_pressure_split.legend(loc='best')
    fig_particle_binder_pressure_split.tight_layout()
    fname = fig_dir + 'particle_binder_pressure_split.svg'
    plt.savefig(fname)

    # =BINDER PARTICLE CONTACTS_P ======================================================================================
    fig_particle_binder_contacts_p, ax_particle_binder_contacts_p = plt.subplots()
    ax_particle_binder_contacts_p.plot(ref.loading_surface_position * 1E2,
                                       ref.loading_particle_contacts_mean,
                                       'k-', zorder=1)
    ax_particle_binder_contacts_p.plot(ref.unloading_surface_position * 1E2,
                                       ref.unloading_particle_contacts_mean,
                                       'k-', zorder=1)
    ax_particle_binder_contacts_p.fill_between(
        ref.loading_surface_position * 1E2,
        (ref.loading_particle_contacts_mean - ref.loading_particle_contacts_std),
        (ref.loading_particle_contacts_mean + ref.loading_particle_contacts_std),
        color='C0', label=ref.label, alpha=1, zorder=1)
    ax_particle_binder_contacts_p.fill_between(
        ref.unloading_surface_position * 1E2,
        (ref.unloading_particle_contacts_mean - ref.unloading_particle_contacts_std),
        (ref.unloading_particle_contacts_mean + ref.unloading_particle_contacts_std),
        color='C0', alpha=1, zorder=1)

    ax_particle_binder_contacts_p.plot(sf_2e8_ef_1e_1.loading_surface_position * 1E2,
                                       sf_2e8_ef_1e_1.loading_particle_contacts_mean,
                                       'k-', zorder=2)
    ax_particle_binder_contacts_p.plot(sf_2e8_ef_1e_1.unloading_surface_position * 1E2,
                                       sf_2e8_ef_1e_1.unloading_particle_contacts_mean,
                                       'k-', zorder=2)
    ax_particle_binder_contacts_p.fill_between(
        sf_2e8_ef_1e_1.loading_surface_position * 1E2,
        (sf_2e8_ef_1e_1.loading_particle_contacts_mean - sf_2e8_ef_1e_1.loading_particle_contacts_std),
        (sf_2e8_ef_1e_1.loading_particle_contacts_mean + sf_2e8_ef_1e_1.loading_particle_contacts_std),
        color='C3', label=sf_2e8_ef_1e_1.label, alpha=1, zorder=2)
    ax_particle_binder_contacts_p.fill_between(
        sf_2e8_ef_1e_1.unloading_surface_position * 1E2,
        (sf_2e8_ef_1e_1.unloading_particle_contacts_mean - sf_2e8_ef_1e_1.unloading_particle_contacts_std),
        (sf_2e8_ef_1e_1.unloading_particle_contacts_mean + sf_2e8_ef_1e_1.unloading_particle_contacts_std),
        color='C3', alpha=1, zorder=2)

    ax_particle_binder_contacts_p.set_xlim(xmin=104.8, xmax=135)
    ax_particle_binder_contacts_p.set_ylim(ymin=0, ymax=6000)
    ax_particle_binder_contacts_p.xaxis.set_major_locator(MultipleLocator(5))
    ax_particle_binder_contacts_p.yaxis.set_major_locator(MultipleLocator(500))
    ax_particle_binder_contacts_p.set_ylabel('Particle contacts [-]')
    ax_particle_binder_contacts_p.set_xlabel('Calendering surface height [µm]')
    ax_particle_binder_contacts_p.legend(loc='best')
    fig_particle_binder_contacts_p.tight_layout()
    fname = fig_dir + 'particle_binder_contacts_p.svg'
    plt.savefig(fname)

    # =BINDER PARICLE CONTACTS_P SPLIT==================================================================================
    fig_particle_binder_contacts_p_split, ax_particle_binder_contacts_p_split = plt.subplots()
    ax_particle_binder_contacts_p_split.set_xlabel('Calendering surface height [µm]')
    ax_particle_binder_contacts_p_split.set_ylabel('Particle contacts [-]')
    # ax_particle_binder_contacts_p_split.set_xlim(xmin=0)
    # ax_particle_binder_contacts_p_split.set_ylim(ymin=0, ymax=180)
    ax_particle_binder_contacts_p_split.plot(ref.loading_surface_position, ref.loading_particle_contacts_matrix.T,
                                             color='C0', label=ref.label)
    ax_particle_binder_contacts_p_split.plot(ref.unloading_surface_position,
                                             ref.unloading_particle_contacts_matrix.T, color='C0')

    ax_particle_binder_contacts_p_split.plot(sf_2e8_ef_1e_1.loading_surface_position,
                                             sf_2e8_ef_1e_1.loading_particle_contacts_matrix.T,
                                             color='C3', label=sf_2e8_ef_1e_1.label)
    ax_particle_binder_contacts_p_split.plot(sf_2e8_ef_1e_1.unloading_surface_position,
                                             sf_2e8_ef_1e_1.unloading_particle_contacts_matrix.T, color='C3')

    ax_particle_binder_contacts_p_split.legend(loc='best')
    fig_particle_binder_contacts_p_split.tight_layout()
    fname = fig_dir + 'particle_binder_contacts_p_split.svg'
    plt.savefig(fname)

    # =BINDER PARTICLE CONTACTS_B=======================================================================================
    fig_particle_binder_contacts_b, ax_particle_binder_contacts_b = plt.subplots()
    ax_particle_binder_contacts_b.plot(ref.loading_surface_position * 1E2,
                                       ref.loading_binder_contacts_mean,
                                       'k-', zorder=1)
    ax_particle_binder_contacts_b.plot(ref.unloading_surface_position * 1E2,
                                       ref.unloading_binder_contacts_mean,
                                       'k-', zorder=1)
    ax_particle_binder_contacts_b.fill_between(
        ref.loading_surface_position * 1E2,
        (ref.loading_binder_contacts_mean - ref.loading_binder_contacts_std),
        (ref.loading_binder_contacts_mean + ref.loading_binder_contacts_std),
        color='C0', label=ref.label, alpha=1, zorder=1)
    ax_particle_binder_contacts_b.fill_between(
        ref.unloading_surface_position * 1E2,
        (ref.unloading_binder_contacts_mean - ref.unloading_binder_contacts_std),
        (ref.unloading_binder_contacts_mean + ref.unloading_binder_contacts_std),
        color='C0', alpha=1, zorder=1)

    ax_particle_binder_contacts_b.plot(sf_2e8_ef_1e_1.loading_surface_position * 1E2,
                                       sf_2e8_ef_1e_1.loading_binder_contacts_mean,
                                       'k-', zorder=2)
    ax_particle_binder_contacts_b.plot(sf_2e8_ef_1e_1.unloading_surface_position * 1E2,
                                       sf_2e8_ef_1e_1.unloading_binder_contacts_mean,
                                       'k-', zorder=2)
    ax_particle_binder_contacts_b.fill_between(
        sf_2e8_ef_1e_1.loading_surface_position * 1E2,
        (sf_2e8_ef_1e_1.loading_binder_contacts_mean - sf_2e8_ef_1e_1.loading_binder_contacts_std),
        (sf_2e8_ef_1e_1.loading_binder_contacts_mean + sf_2e8_ef_1e_1.loading_binder_contacts_std),
        color='C3', label=sf_2e8_ef_1e_1.label, alpha=1, zorder=2)
    ax_particle_binder_contacts_b.fill_between(
        sf_2e8_ef_1e_1.unloading_surface_position * 1E2,
        (sf_2e8_ef_1e_1.unloading_binder_contacts_mean - sf_2e8_ef_1e_1.unloading_binder_contacts_std),
        (sf_2e8_ef_1e_1.unloading_binder_contacts_mean + sf_2e8_ef_1e_1.unloading_binder_contacts_std),
        color='C3', alpha=1, zorder=2)

    ax_particle_binder_contacts_b.set_xlim(xmin=104.8, xmax=135)
    ax_particle_binder_contacts_b.set_ylim(ymin=8800, ymax=10200)
    ax_particle_binder_contacts_b.xaxis.set_major_locator(MultipleLocator(5))
    ax_particle_binder_contacts_b.yaxis.set_major_locator(MultipleLocator(200))
    ax_particle_binder_contacts_b.set_ylabel('Binder contacts [-]')
    ax_particle_binder_contacts_b.set_xlabel('Calendering surface height [µm]')
    fig_particle_binder_contacts_b.tight_layout()
    ax_particle_binder_contacts_b.legend(loc='best')
    fname = fig_dir + 'particle_binder_contacts_b.svg'
    plt.savefig(fname)

    # =BINDER PARTICLE CONTACTS_B SPLIT=================================================================================
    fig_particle_binder_contacts_b_split, ax_particle_binder_contacts_b_split = plt.subplots()
    ax_particle_binder_contacts_b_split.set_xlabel('Calendering surface height [µm]')
    ax_particle_binder_contacts_b_split.set_ylabel('Binder contacts [-]')
    # ax_particle_binder_contacts_b_split.set_xlim(xmin=0)
    # ax_particle_binder_contacts_b_split.set_ylim(ymin=0, ymax=180)
    ax_particle_binder_contacts_b_split.plot(ref.loading_surface_position, ref.loading_binder_contacts_matrix.T,
                                             color='C0', label=ref.label)
    ax_particle_binder_contacts_b_split.plot(ref.unloading_surface_position,
                                             ref.unloading_binder_contacts_matrix.T, color='C0')

    ax_particle_binder_contacts_b_split.plot(sf_2e8_ef_1e_1.loading_surface_position,
                                             sf_2e8_ef_1e_1.loading_binder_contacts_matrix.T,
                                             color='C3', label=sf_2e8_ef_1e_1.label)
    ax_particle_binder_contacts_b_split.plot(sf_2e8_ef_1e_1.unloading_surface_position,
                                             sf_2e8_ef_1e_1.unloading_binder_contacts_matrix.T, color='C3')

    ax_particle_binder_contacts_b_split.legend(loc='best')
    fig_particle_binder_contacts_b_split.tight_layout()
    fname = fig_dir + 'particle_binder_contacts_b_split.svg'
    plt.savefig(fname)

    # =BINDER PARTICLE FRACTURE_P ======================================================================================
    fig_particle_binder_fracture_p, ax_particle_binder_fracture_p = plt.subplots()

    ax_particle_binder_fracture_p.plot(sf_2e8_ef_1e_1.loading_surface_position * 1E2,
                                       sf_2e8_ef_1e_1.loading_fractured_particles_mean * 1E2,
                                       'k-', zorder=1)
    ax_particle_binder_fracture_p.plot(sf_2e8_ef_1e_1.unloading_surface_position * 1E2,
                                       sf_2e8_ef_1e_1.unloading_fractured_particles_mean * 1E2,
                                       'k-', zorder=1)
    ax_particle_binder_fracture_p.fill_between(
        sf_2e8_ef_1e_1.loading_surface_position * 1E2,
        (sf_2e8_ef_1e_1.loading_fractured_particles_mean - sf_2e8_ef_1e_1.loading_fractured_particles_std) * 1E2,
        (sf_2e8_ef_1e_1.loading_fractured_particles_mean + sf_2e8_ef_1e_1.loading_fractured_particles_std) * 1E2,
        color='C3', label=sf_2e8_ef_1e_1.label, alpha=1, zorder=1)
    ax_particle_binder_fracture_p.fill_between(
        sf_2e8_ef_1e_1.unloading_surface_position * 1E2,
        (sf_2e8_ef_1e_1.unloading_fractured_particles_mean - sf_2e8_ef_1e_1.unloading_fractured_particles_std) * 1E2,
        (sf_2e8_ef_1e_1.unloading_fractured_particles_mean + sf_2e8_ef_1e_1.unloading_fractured_particles_std) * 1E2,
        color='C3', alpha=1, zorder=1)

    ax_particle_binder_fracture_p.set_xlim(xmin=104.8, xmax=135)
    ax_particle_binder_fracture_p.xaxis.set_major_locator(MultipleLocator(5))
    ax_particle_binder_fracture_p.set_ylim(ymin=0, ymax=50)
    ax_particle_binder_fracture_p.yaxis.set_major_locator(MultipleLocator(10))
    ax_particle_binder_fracture_p.set_ylabel('Fractured particles [%]')
    ax_particle_binder_fracture_p.set_xlabel('Calendering surface height [µm]')
    ax_particle_binder_fracture_p.legend(loc='best')
    fig_particle_binder_fracture_p.tight_layout()
    fname = fig_dir + 'particle_binder_fracture_p.svg'
    plt.savefig(fname)

    # =BINDER PARTICLE FRACTURE_P SPLIT=================================================================================
    fig_particle_binder_fracture_p_split, ax_particle_binder_fracture_p_split = plt.subplots()
    ax_particle_binder_fracture_p_split.set_xlabel('Calendering surface height [µm]')
    ax_particle_binder_fracture_p_split.set_ylabel('Fraction of fractured particles [-]')
    # ax_particle_binder_fracture_p_split.set_xlim(xmin=0)
    # ax_particle_binder_contacts_p_split.set_ylim(ymin=0, ymax=180)
    ax_particle_binder_fracture_p_split.plot(sf_2e8_ef_1e_1.loading_surface_position,
                                             sf_2e8_ef_1e_1.loading_fractured_particles_matrix.T,
                                             color='C3', label=sf_2e8_ef_1e_1.label)
    ax_particle_binder_fracture_p_split.plot(sf_2e8_ef_1e_1.unloading_surface_position,
                                             sf_2e8_ef_1e_1.unloading_fractured_particles_matrix.T, color='C3')

    ax_particle_binder_fracture_p_split.legend(loc='best')
    fig_particle_binder_fracture_p_split.tight_layout()
    fname = fig_dir + 'particle_binder_fracture_p_split.svg'
    plt.savefig(fname)

    # =BINDER PARTICLE FRACTURE_B=======================================================================================
    fig_particle_binder_fracture_b, ax_particle_binder_fracture_b = plt.subplots()

    ax_particle_binder_fracture_b.plot(sf_2e8_ef_1e_1.loading_surface_position * 1E2,
                                       sf_2e8_ef_1e_1.loading_fractured_binders_mean * 1E2,
                                       'k-', zorder=1)
    ax_particle_binder_fracture_b.plot(sf_2e8_ef_1e_1.unloading_surface_position * 1E2,
                                       sf_2e8_ef_1e_1.unloading_fractured_binders_mean * 1E2,
                                       'k-', zorder=1)
    ax_particle_binder_fracture_b.fill_between(
        sf_2e8_ef_1e_1.loading_surface_position * 1E2,
        (sf_2e8_ef_1e_1.loading_fractured_binders_mean - sf_2e8_ef_1e_1.loading_fractured_binders_std) * 1E2,
        (sf_2e8_ef_1e_1.loading_fractured_binders_mean + sf_2e8_ef_1e_1.loading_fractured_binders_std) * 1E2,
        color='C3', label=sf_2e8_ef_1e_1.label, alpha=1, zorder=1)
    ax_particle_binder_fracture_b.fill_between(
        sf_2e8_ef_1e_1.unloading_surface_position * 1E2,
        (sf_2e8_ef_1e_1.unloading_fractured_binders_mean - sf_2e8_ef_1e_1.unloading_fractured_binders_std) * 1E2,
        (sf_2e8_ef_1e_1.unloading_fractured_binders_mean + sf_2e8_ef_1e_1.unloading_fractured_binders_std) * 1E2,
        color='C3', alpha=1, zorder=1)

    ax_particle_binder_fracture_b.set_xlim(xmin=104.8, xmax=135)
    ax_particle_binder_fracture_b.xaxis.set_major_locator(MultipleLocator(5))
    ax_particle_binder_fracture_b.set_ylim(ymin=0, ymax=5)
    ax_particle_binder_fracture_b.yaxis.set_major_locator(MultipleLocator(1))
    ax_particle_binder_fracture_b.set_ylabel('Fractured binder contacts [%]')
    ax_particle_binder_fracture_b.set_xlabel('Calendering surface height [µm]')
    fig_particle_binder_fracture_b.tight_layout()
    ax_particle_binder_fracture_b.legend(loc='best')
    fname = fig_dir + 'particle_binder_fracture_b.svg'
    plt.savefig(fname)

    # =BINDER PARTICLE FRACTURE_B SPLIT=================================================================================
    fig_particle_binder_fracture_b_split, ax_particle_binder_fracture_b_split = plt.subplots()
    ax_particle_binder_fracture_b_split.set_xlabel('Calendering surface height [µm]')
    ax_particle_binder_fracture_b_split.set_ylabel('Fraction of fractured binder contacts [-]')
    # ax_particle_binder_fracture_b_split.set_xlim(xmin=0)
    # ax_particle_binder_fracture_b_split.set_ylim(ymin=0, ymax=180)

    ax_particle_binder_fracture_b_split.plot(sf_2e8_ef_1e_1.loading_surface_position,
                                             sf_2e8_ef_1e_1.loading_fractured_binders_matrix.T,
                                             color='C3', label=sf_2e8_ef_1e_1.label)
    ax_particle_binder_fracture_b_split.plot(sf_2e8_ef_1e_1.unloading_surface_position,
                                             sf_2e8_ef_1e_1.unloading_fractured_binders_matrix.T, color='C3')

    fig_particle_binder_fracture_b_split.tight_layout()
    fname = fig_dir + 'particle_binder_fracture_b_split.svg'
    plt.savefig(fname)

    ########################################################################################################################
    plt.show()
