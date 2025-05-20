#
# Created by Axel on 2025-05-20
# Adapted from article_4/cycling_spreads.py and article_3/in_plane_stress_spreads.py
#
from Bertil_functions.Bertil_functions import *
from force_model_impact_on_calendering.Bertil_calendering_pressure import local_data_gatherer, bertil_data_gatherer, \
    contact_counter_bertil
from force_model_impact_on_calendering.Bertil_calendering_pressure_multiple_simulations import \
    calendering_plot_processing
from force_model_impact_on_calendering.Bertil_fractured_particles import fractured_particle_gatherer
from force_model_impact_on_calendering.Bertil_fractured_binder_contacts import fractured_binder_counter
from force_model_impact_on_calendering.Bertil_mechanical_properties_multiple_runs import mech_plot_prop, stiffness_func
from force_model_impact_on_calendering.Bertil_mechanical_properties import stress_and_linear_strain_finder

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import os
import matplotlib
import shutil
from scipy import interpolate

matplotlib.style.use('axel_style')


def rve_stress_calculator(periodic_BC_data, force_fabric_tensor_data, layer_height_data):
    time = periodic_BC_data[:, 0]

    x_side_length = periodic_BC_data[:, 2] - periodic_BC_data[:, 1]
    y_side_length = periodic_BC_data[:, 4] - periodic_BC_data[:, 3]
    print('x and y side length is = ' + str(x_side_length[0]) + ', ' + str(y_side_length[0]))

    print('layer height = ' + str(layer_height_data[0]))
    vol = x_side_length * y_side_length * layer_height_data
    sxx = -force_fabric_tensor_data[:, 1] / vol
    syy = -force_fabric_tensor_data[:, 5] / vol
    szz = -force_fabric_tensor_data[:, 9] / vol

    tau_xy = -force_fabric_tensor_data[:, 2] / vol
    tau_xz = -force_fabric_tensor_data[:, 3] / vol

    tau_yx = -force_fabric_tensor_data[:, 4] / vol
    tau_yz = -force_fabric_tensor_data[:, 6] / vol

    tau_zx = -force_fabric_tensor_data[:, 7] / vol
    tau_zy = -force_fabric_tensor_data[:, 8] / vol

    exx = (x_side_length - x_side_length[0])/x_side_length[0]

    return time, vol, exx, sxx, syy, szz, tau_xy, tau_xz, tau_yx, tau_yz, tau_zx, tau_zy


def in_plane_stress_processing(sim_case, no_sims):
    simulation_dictionary = {}
    for i in range(1, no_sims + 1):
        sim_dir = '/scratch/users/axlun/DEMsim/results/article_4/final_runs/' + str(i) + '/' + sim_case + \
                  '/electrode_uniaxial_tension'
        print(sim_dir)
        simulation_dictionary[i] = Cycling(sim_dir)
    time_size = simulation_dictionary[1].time.size
    sxx_matrix = np.zeros((no_sims, time_size))
    syy_matrix = np.zeros((no_sims, time_size))
    sxy_matrix = np.zeros((no_sims, time_size))
    pressure_matrix = np.zeros((no_sims, time_size))
    height_matrix = np.zeros((no_sims, time_size))
    binder_contacts_matrix = np.zeros((no_sims, time_size))
    initialised_binder_contacts_matrix = np.zeros((no_sims, time_size))
    fractured_binders_matrix = np.zeros((no_sims, time_size))
    particle_contacts_matrix = np.zeros((no_sims, time_size))
    fractured_particles_matrix = np.zeros((no_sims, time_size))
    time = simulation_dictionary[1].absolute_time
    exx = simulation_dictionary[1].exx

    for i in range(no_sims):
        sxx_matrix[i, :] = simulation_dictionary[i + 1].sxx
        syy_matrix[i, :] = simulation_dictionary[i + 1].syy
        sxy_matrix[i, :] = simulation_dictionary[i + 1].tau_xy
        pressure_matrix[i, :] = simulation_dictionary[i + 1].avg_in_plane_compression
        binder_contacts_matrix[i, :] = simulation_dictionary[i + 1].active_binder_contact_vec
        initialised_binder_contacts_matrix[i, :] = simulation_dictionary[i + 1].initialised_binder_contact_vec
        fractured_binders_matrix[i, :] = simulation_dictionary[i + 1].fractured_binder_contact_vec
        particle_contacts_matrix[i, :] = simulation_dictionary[i + 1].particle_contact_vec
        fractured_particles_matrix[i, :] = simulation_dictionary[i + 1].fractured_particles

    # =NORMALISATION OF FRACTURE MATRIX=================================================================================
    fractured_binders_matrix = fractured_binders_matrix / initialised_binder_contacts_matrix
    # fractured_binders_matrix = fractured_binders_matrix / binder_contacts_matrix
    fractured_particles_matrix = fractured_particles_matrix / 5000
    # ==================================================================================================================
    sxx_mean = np.mean(sxx_matrix, axis=0)
    sxx_std = np.std(sxx_matrix, axis=0)

    syy_mean = np.mean(syy_matrix, axis=0)
    syy_std = np.std(syy_matrix, axis=0)

    sxy_mean = np.mean(sxy_matrix, axis=0)
    sxy_std = np.std(sxy_matrix, axis=0)

    in_plane_stress_mean = np.mean(pressure_matrix, axis=0)
    in_plane_stress_std = np.std(pressure_matrix, axis=0)
    in_plane_stress_matrix = pressure_matrix

    binder_contacts_mean = np.mean(binder_contacts_matrix, axis=0)
    binder_contacts_std = np.std(binder_contacts_matrix, axis=0)

    fractured_binders_mean = np.mean(fractured_binders_matrix, axis=0)
    fractured_binders_std = np.std(fractured_binders_matrix, axis=0)

    particle_contacts_mean = np.mean(particle_contacts_matrix, axis=0)
    particle_contacts_std = np.std(particle_contacts_matrix, axis=0)

    fractured_particles_mean = np.mean(fractured_particles_matrix, axis=0)
    fractured_particles_std = np.std(fractured_particles_matrix, axis=0)

    height_mean = np.mean(height_matrix, axis=0)
    height_std = np.std(height_matrix, axis=0)

    return time, exx, sxx_mean, sxx_std, sxx_matrix, \
           syy_mean, syy_std, syy_matrix, \
           sxy_mean, sxy_std, sxy_matrix, \
           in_plane_stress_mean, in_plane_stress_std, in_plane_stress_matrix, \
           particle_contacts_mean, particle_contacts_std, particle_contacts_matrix, \
           binder_contacts_mean, binder_contacts_std, binder_contacts_matrix, \
           fractured_particles_mean, fractured_particles_std, fractured_particles_matrix, \
           fractured_binders_mean, fractured_binders_std, fractured_binders_matrix, \
           height_mean, height_std, height_matrix


class Cycling:
    def __init__(self, sim_dir):
        self.sim_dir = sim_dir
        self.force_data, self.surface_force_index, self.surface_position_index, \
        self.surface_position_data, self.periodic_BC_data, self.force_fabric_tensor_data, \
        self.ke = bertil_data_gatherer(sim_dir)

        self.height_time, self.avg_height = bertil_layer_height(sim_dir, 50)

        self.time, self.rve_vol, self.exx, self.sxx, self.syy, self.szz, self.tau_xy, self.tau_xz, self.tau_yx, \
        self.tau_yz, self.tau_zx, self.tau_zy = \
            rve_stress_calculator(self.periodic_BC_data, self.force_fabric_tensor_data, self.avg_height)

        self.absolute_time = self.time - self.time[0]
        self.avg_in_plane_compression = - (self.sxx + self.syy) / 2

        # Also get the number of contacts here
        self.contacts_time, self.particle_contact_vec, self.binder_contact_vec, self.binder_particle_contact_vec \
            = contact_counter_bertil(sim_dir)

        # As well as the number of fractures
        self.frac_binder_time, self.initialised_binder_contact_vec, self.active_binder_contact_vec, \
        self.fractured_binder_contact_vec = fractured_binder_counter(sim_dir)

        self.frac_particle_time, self.fractured_particles, self.fracture_particle_array = \
            fractured_particle_gatherer(sim_dir)
        # ===============================================================================================================


class Simulation:
    def __init__(self, sim_case, no_sims, label='template'):
        self.time, self.exx, self.sxx_mean, self.sxx_std, self.sxx_matrix, \
        self.syy_mean, self.syy_std, self.syy_matrix, \
        self.sxy_mean, self.sxy_std, self.sxy_matrix, \
        self.in_plane_stress_mean, self.in_plane_stress_std, self.in_plane_stress_matrix, \
        self.particle_contacts_mean, self.particle_contacts_std, self.particle_contacts_matrix, \
        self.binder_contacts_mean, self.binder_contacts_std, self.binder_contacts_matrix, \
        self.fractured_particles_mean, self.fractured_particles_std, self.fractured_particles_matrix, \
        self.fractured_binders_mean, self.fractured_binders_std, self.fractured_binders_matrix, \
        self.height_mean, self.height_std, self.height_matrix \
            = in_plane_stress_processing(sim_case, no_sims)

        self.label = label


if __name__ == '__main__':
    # ==================================================================================================================
    # TODO:
    #       *
    # ==================================================================================================================

    # ==FIGURE SAVE DIR=================================================================================================
    fig_dir = 'C:/temp/figures/article_4/final_runs/tensile_spreads/'
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

    rest_time = 3
    end_time = ref.time[-1]- rest_time
    end_strain = ref.exx[-1]
    no_strain_ticks = 15
    tick_array = np.linspace(0,end_time, no_strain_ticks)
    tick_labels = np.linspace(0, end_strain, no_strain_ticks)
    end_plot_strain = 15
    # n_cycles = 3
    # n_soc_ticks = n_cycles * 2 + 1
    # tick_array = np.linspace(0,end_plot_strain,n_soc_ticks)
    # soc_tick_labels = [0, 100, 0, 100, 0, 100, 0]
    ########################################PARTICLE FRACTURE###########################################################
    # =FIGURE PARTICLE PRESSURE=========================================================================================
    fig_particle_pressure, ax_particle_pressure = plt.subplots()
    ax_particle_pressure.plot(ref.exx * 1E2,
                              ref.sxx_mean* 1E-6,
                              'k-', zorder=1)
    ax_particle_pressure.fill_between(
        ref.exx * 1E2,
        (ref.sxx_mean - ref.sxx_std) * 1E-6,
        (ref.sxx_mean + ref.sxx_std) * 1E-6,
        color='C0', label=ref.label, alpha=1, zorder=1)

    ax_particle_pressure.plot(sf_2e8.exx * 1E2,
                              sf_2e8.sxx_mean * 1E-6,
                              'k-', zorder=2)
    ax_particle_pressure.fill_between(
        sf_2e8.exx * 1E2,
        (sf_2e8.sxx_mean - sf_2e8.sxx_std) * 1E-6,
        (sf_2e8.sxx_mean + sf_2e8.sxx_std) * 1E-6,
        color='C1', label=sf_2e8.label, alpha=1, zorder=2)

    ax_particle_pressure.plot(sf_1e8.exx * 1E2,
                              sf_1e8.sxx_mean * 1E-6,
                              'k-', zorder=3)
    ax_particle_pressure.fill_between(
        sf_1e8.exx * 1E2,
        (sf_1e8.sxx_mean - sf_1e8.sxx_std) * 1E-6,
        (sf_1e8.sxx_mean + sf_1e8.sxx_std) * 1E-6,
        color='C2', label=sf_1e8.label, alpha=1, zorder=3)

    ax_particle_pressure.plot(sf_5e7.exx * 1E2,
                              sf_5e7.sxx_mean * 1E-6,
                              'k-', zorder=4)
    ax_particle_pressure.fill_between(
        sf_5e7.exx * 1E2,
        (sf_5e7.sxx_mean - sf_5e7.sxx_std) * 1E-6,
        (sf_5e7.sxx_mean + sf_5e7.sxx_std) * 1E-6,
        color='C3', label=sf_5e7.label, alpha=1, zorder=4)
    ax_particle_pressure.set_xlim(xmin=0, xmax=end_plot_strain)
    ax_particle_pressure.set_ylim(ymin=-20, ymax=120)
    ax_particle_pressure.yaxis.set_major_locator(MultipleLocator(20))
    ax_particle_pressure.xaxis.set_major_locator(MultipleLocator(1))
    ax_particle_pressure.set_ylabel('Stress [MPa]')
    ax_particle_pressure.set_xlabel('Strain [%]')

    fig_particle_pressure.tight_layout()
    ax_particle_pressure.legend(loc='best')
    fname = fig_dir + 'particle_pressure.svg'
    plt.savefig(fname)

    # =PARTICLE PRESSURE SPLIT==========================================================================================
    fig_particle_pressure_split, ax_particle_pressure_split = plt.subplots()
    ax_particle_pressure_split.set_xlabel('Strain [%]')
    ax_particle_pressure_split.set_ylabel('Stress [MPa]')
    # ax_particle_pressure_split.set_xlim(xmin=0)
    # ax_particle_pressure_split.set_ylim(ymin=0, ymax=180)


    ax_particle_pressure_split.plot(ref.exx * 1E2, ref.sxx_matrix.T * 1E-6,
                                    color='C0', label=ref.label)
    ax_particle_pressure_split.plot(ref.exx * 1E2, ref.sxx_matrix.T * 1E-6,
                                    color='C0', label=ref.label)

    ax_particle_pressure_split.plot(sf_5e7.exx * 1E2, sf_5e7.sxx_matrix.T * 1E-6,
                                    color='C1', label=sf_5e7.label)

    ax_particle_pressure_split.plot(sf_1e8.exx * 1E2, sf_1e8.sxx_matrix.T * 1E-6,
                                    color='C2', label=sf_1e8.label)

    ax_particle_pressure_split.plot(sf_2e8.exx * 1E2, sf_2e8.sxx_matrix.T * 1E-6,
                                    color='C3', label=sf_2e8.label)
    ax_particle_pressure_split.legend(loc='best')
    fig_particle_pressure_split.tight_layout()
    fname = fig_dir + 'particle_pressure_split.svg'
    plt.savefig(fname)

    # =PARTICLE CONTACTS_P ===============================================================================================
    fig_particle_contacts_p, ax_particle_contacts_p = plt.subplots()
    ax_particle_contacts_p.plot(ref.exx * 1E2,
                                ref.particle_contacts_mean,
                                'k-', zorder=1)
    ax_particle_contacts_p.fill_between(
        ref.exx * 1E2,
        (ref.particle_contacts_mean - ref.particle_contacts_std),
        (ref.particle_contacts_mean + ref.particle_contacts_std),
        color='C0', label=ref.label, alpha=1, zorder=1)

    ax_particle_contacts_p.plot(sf_2e8.exx * 1E2,
                                sf_2e8.particle_contacts_mean,
                                'k-', zorder=2)
    ax_particle_contacts_p.fill_between(
        sf_2e8.exx * 1E2,
        (sf_2e8.particle_contacts_mean - sf_2e8.particle_contacts_std),
        (sf_2e8.particle_contacts_mean + sf_2e8.particle_contacts_std),
        color='C1', label=sf_2e8.label, alpha=1, zorder=2)

    ax_particle_contacts_p.plot(sf_1e8.exx * 1E2,
                                sf_1e8.particle_contacts_mean,
                                'k-', zorder=3)
    ax_particle_contacts_p.fill_between(
        sf_1e8.exx * 1E2,
        (sf_1e8.particle_contacts_mean - sf_1e8.particle_contacts_std),
        (sf_1e8.particle_contacts_mean + sf_1e8.particle_contacts_std),
        color='C2', label=sf_1e8.label, alpha=1, zorder=3)

    ax_particle_contacts_p.plot(sf_5e7.exx * 1E2,
                                sf_5e7.particle_contacts_mean,
                                'k-', zorder=4)
    ax_particle_contacts_p.fill_between(
        sf_5e7.exx * 1E2,
        (sf_5e7.particle_contacts_mean - sf_5e7.particle_contacts_std),
        (sf_5e7.particle_contacts_mean + sf_5e7.particle_contacts_std),
        color='C3', label=sf_5e7.label, alpha=1, zorder=4)

    ax_particle_contacts_p.set_xlim(xmin=0, xmax=end_plot_strain)
    ax_particle_contacts_p.set_ylim(ymin=2000, ymax=4000)
    ax_particle_contacts_p.yaxis.set_major_locator(MultipleLocator(250))
    ax_particle_contacts_p.xaxis.set_major_locator(MultipleLocator(1))
    ax_particle_contacts_p.set_ylabel('Particle contacts [-]')
    ax_particle_contacts_p.set_xlabel('Strain [%]')
    ax_particle_contacts_p.legend(loc='best')

    fig_particle_contacts_p.tight_layout()
    fname = fig_dir + 'particle_contacts_p.svg'
    plt.savefig(fname)

    # =PARTICLE CONTACTS_P SPLIT========================================================================================
    fig_particle_contacts_p_split, ax_particle_contacts_p_split = plt.subplots()
    ax_particle_contacts_p_split.set_xlabel('Strain [%]')
    ax_particle_contacts_p_split.set_ylabel('Particle contacts [-]')
    # ax_particle_contacts_p_split.set_xlim(xmin=0)
    # ax_particle_contacts_p_split.set_ylim(ymin=0, ymax=180)
    ax_particle_contacts_p_split.plot(ref.exx * 1E2, ref.particle_contacts_matrix.T,
                                      color='C0', label=ref.label)

    ax_particle_contacts_p_split.plot(sf_5e7.exx * 1E2, sf_5e7.particle_contacts_matrix.T,
                                      color='C1', label=sf_5e7.label)

    ax_particle_contacts_p_split.plot(sf_1e8.exx * 1E2, sf_1e8.particle_contacts_matrix.T,
                                      color='C2', label=sf_1e8.label)

    ax_particle_contacts_p_split.plot(sf_2e8.exx * 1E2, sf_2e8.particle_contacts_matrix.T,
                                      color='C3', label=sf_2e8.label)
    ax_particle_contacts_p_split.legend(loc='best')
    fig_particle_contacts_p_split.tight_layout()
    fname = fig_dir + 'particle_contacts_p_split.svg'
    plt.savefig(fname)

    # =PARTICLE CONTACTS_B==============================================================================================
    fig_particle_contacts_b, ax_particle_contacts_b = plt.subplots()
    ax_particle_contacts_b.plot(ref.exx * 1E2,
                                ref.binder_contacts_mean,
                                'k-', zorder=1)
    ax_particle_contacts_b.fill_between(
        ref.exx * 1E2,
        (ref.binder_contacts_mean - ref.binder_contacts_std),
        (ref.binder_contacts_mean + ref.binder_contacts_std),
        color='C0', label=ref.label, alpha=1, zorder=1)

    ax_particle_contacts_b.plot(sf_2e8.exx * 1E2,
                                sf_2e8.binder_contacts_mean,
                                'k-', zorder=2)
    ax_particle_contacts_b.fill_between(
        sf_2e8.exx * 1E2,
        (sf_2e8.binder_contacts_mean - sf_2e8.binder_contacts_std),
        (sf_2e8.binder_contacts_mean + sf_2e8.binder_contacts_std),
        color='C1', label=sf_2e8.label, alpha=1, zorder=2)

    ax_particle_contacts_b.plot(sf_1e8.exx * 1E2,
                                sf_1e8.binder_contacts_mean,
                                'k-', zorder=3)
    ax_particle_contacts_b.fill_between(
        sf_1e8.exx * 1E2,
        (sf_1e8.binder_contacts_mean - sf_1e8.binder_contacts_std),
        (sf_1e8.binder_contacts_mean + sf_1e8.binder_contacts_std),
        color='C2', label=sf_1e8.label, alpha=1, zorder=3)

    ax_particle_contacts_b.plot(sf_5e7.exx * 1E2,
                                sf_5e7.binder_contacts_mean,
                                'k-', zorder=4)
    ax_particle_contacts_b.fill_between(
        sf_5e7.exx * 1E2,
        (sf_5e7.binder_contacts_mean - sf_5e7.binder_contacts_std),
        (sf_5e7.binder_contacts_mean + sf_5e7.binder_contacts_std),
        color='C3', label=sf_5e7.label, alpha=1, zorder=4)

    ax_particle_contacts_b.set_xlim(xmin=0, xmax=end_plot_strain)
    ax_particle_contacts_b.set_ylim(ymin=9500, ymax=10000)
    ax_particle_contacts_b.yaxis.set_major_locator(MultipleLocator(50))
    ax_particle_contacts_b.xaxis.set_major_locator(MultipleLocator(1))
    ax_particle_contacts_b.set_ylabel('Binder contacts [-]')
    ax_particle_contacts_b.set_xlabel('Strain [%]')


    fig_particle_contacts_b.tight_layout()
    ax_particle_contacts_b.legend(loc='best')
    fname = fig_dir + 'particle_contacts_b.svg'
    plt.savefig(fname)

    # =PARTICLE CONTACTS_B SPLIT========================================================================================
    fig_particle_contacts_b_split, ax_particle_contacts_b_split = plt.subplots()
    ax_particle_contacts_b_split.set_xlabel('Strain [%]')
    ax_particle_contacts_b_split.set_ylabel('Binder contacts [-]')
    ax_particle_contacts_b_split.plot(ref.exx * 1E2, ref.binder_contacts_matrix.T,
                                      color='C0', label=ref.label)

    ax_particle_contacts_b_split.plot(sf_5e7.exx * 1E2, sf_5e7.binder_contacts_matrix.T,
                                      color='C1', label=sf_5e7.label)
    ax_particle_contacts_b_split.plot(sf_1e8.exx * 1E2, sf_1e8.binder_contacts_matrix.T,
                                      color='C2', label=sf_1e8.label)

    ax_particle_contacts_b_split.plot(sf_2e8.exx * 1E2, sf_2e8.binder_contacts_matrix.T,
                                      color='C3', label=sf_2e8.label)
    ax_particle_contacts_b_split.legend(loc='best')
    fig_particle_contacts_b_split.tight_layout()
    fname = fig_dir + 'particle_contacts_b_split.svg'
    plt.savefig(fname)

    # =PARTICLE FRACTURE_P ===============================================================================================
    fig_particle_fracture_p, ax_particle_fracture_p = plt.subplots()
    ax_particle_fracture_p.plot(sf_2e8.exx * 1E2,
                                sf_2e8.fractured_particles_mean* 1E2,
                                'k-', zorder=1)
    ax_particle_fracture_p.fill_between(
        sf_2e8.exx * 1E2,
        (sf_2e8.fractured_particles_mean - sf_2e8.fractured_particles_std)* 1E2,
        (sf_2e8.fractured_particles_mean + sf_2e8.fractured_particles_std)* 1E2,
        color='C1', label=sf_2e8.label, alpha=1, zorder=1)

    ax_particle_fracture_p.plot(sf_1e8.exx * 1E2,
                                sf_1e8.fractured_particles_mean* 1E2,
                                'k-', zorder=2)
    ax_particle_fracture_p.fill_between(
        sf_1e8.exx * 1E2,
        (sf_1e8.fractured_particles_mean - sf_1e8.fractured_particles_std)* 1E2,
        (sf_1e8.fractured_particles_mean + sf_1e8.fractured_particles_std)* 1E2,
        color='C2', label=sf_1e8.label, alpha=1, zorder=2)

    ax_particle_fracture_p.plot(sf_5e7.exx * 1E2,
                                sf_5e7.fractured_particles_mean* 1E2,
                                'k-', zorder=3)
    ax_particle_fracture_p.fill_between(
        sf_5e7.exx * 1E2,
        (sf_5e7.fractured_particles_mean - sf_5e7.fractured_particles_std)* 1E2,
        (sf_5e7.fractured_particles_mean + sf_5e7.fractured_particles_std)* 1E2,
        color='C3', label=sf_5e7.label, alpha=1, zorder=3)

    ax_particle_fracture_p.set_xlim(xmin=0, xmax=end_plot_strain)
    ax_particle_fracture_p.set_ylim(ymin=0, ymax=1E2)
    ax_particle_fracture_p.yaxis.set_major_locator(MultipleLocator(10))
    ax_particle_fracture_p.xaxis.set_major_locator(MultipleLocator(1))
    ax_particle_fracture_p.set_ylabel('Fractured particles [%]')
    ax_particle_fracture_p.set_xlabel('Strain [%]')
    ax_particle_fracture_p.legend(loc='best')


    fig_particle_fracture_p.tight_layout()
    fname = fig_dir + 'particle_fracture_p.svg'
    plt.savefig(fname)

    # =PARTICLE FRACTURE_P SPLIT========================================================================================
    fig_particle_fracture_p_split, ax_particle_fracture_p_split = plt.subplots()
    ax_particle_fracture_p_split.set_xlabel('Strain [%]')
    ax_particle_fracture_p_split.set_ylabel('Fraction of fractured particles [-]')
    # ax_particle_fracture_p_split.set_xlim(xmin=0)
    # ax_particle_contacts_p_split.set_ylim(ymin=0, ymax=180)
    ax_particle_fracture_p_split.plot(sf_5e7.exx * 1E2, sf_5e7.fractured_particles_matrix.T,
                                      color='C1', label=sf_5e7.label)
    ax_particle_fracture_p_split.plot(sf_1e8.exx * 1E2, sf_1e8.fractured_particles_matrix.T,
                                      color='C2', label=sf_1e8.label)

    ax_particle_fracture_p_split.plot(sf_2e8.exx * 1E2, sf_2e8.fractured_particles_matrix.T,
                                      color='C3', label=sf_2e8.label)

    ax_particle_fracture_p_split.legend(loc='best')
    fig_particle_fracture_p_split.tight_layout()
    fname = fig_dir + 'particle_fracture_p_split.svg'
    plt.savefig(fname)

    # =PARTICLE FRACTURE_B==============================================================================================
    fig_particle_fracture_b, ax_particle_fracture_b = plt.subplots()
    ax_particle_fracture_b.plot(sf_2e8.exx * 1E2,
                                sf_2e8.fractured_binders_mean* 1E2,
                                'k-', zorder=1)
    ax_particle_fracture_b.fill_between(sf_2e8.exx * 1E2,
                                        (sf_2e8.fractured_binders_mean - sf_2e8.fractured_binders_std)* 1E2,
                                        (sf_2e8.fractured_binders_mean + sf_2e8.fractured_binders_std)* 1E2,
                                        color='C1', label=sf_2e8.label, alpha=1, zorder=1)

    ax_particle_fracture_b.plot(sf_1e8.exx * 1E2,
                                sf_1e8.fractured_binders_mean* 1E2,
                                'k-', zorder=2)
    ax_particle_fracture_b.fill_between(sf_1e8.exx * 1E2,
                                        (sf_1e8.fractured_binders_mean - sf_1e8.fractured_binders_std)* 1E2,
                                        (sf_1e8.fractured_binders_mean + sf_1e8.fractured_binders_std)* 1E2,
                                        color='C2', label=sf_1e8.label, alpha=1, zorder=2)

    ax_particle_fracture_b.plot(sf_5e7.exx * 1E2,
                                sf_5e7.fractured_binders_mean* 1E2,
                                'k-', zorder=3)
    ax_particle_fracture_b.fill_between(sf_5e7.exx * 1E2,
                                        (sf_5e7.fractured_binders_mean - sf_5e7.fractured_binders_std)* 1E2,
                                        (sf_5e7.fractured_binders_mean + sf_5e7.fractured_binders_std)* 1E2,
                                        color='C3', label=sf_5e7.label, alpha=1, zorder=3)

    ax_particle_fracture_b.set_xlim(xmin=0, xmax=end_plot_strain)
    ax_particle_fracture_b.set_ylim(ymin=0)
    ax_particle_fracture_b.set_ylabel('Fractured binder contacts [%]')
    ax_particle_fracture_b.set_xlabel('Strain [%]')


    fig_particle_fracture_b.tight_layout()
    ax_particle_fracture_b.legend(loc='best')
    fname = fig_dir + 'particle_fracture_b.svg'
    plt.savefig(fname)

    # =PARTICLE FRACTURE_B SPLIT========================================================================================
    fig_particle_fracture_b_split, ax_particle_fracture_b_split = plt.subplots()
    ax_particle_fracture_b_split.set_xlabel('Strain [%]')
    ax_particle_fracture_b_split.set_ylabel('Fraction of fractured binder contacts [-]')
    # ax_particle_fracture_b_split.set_xlim(xmin=0)
    # ax_particle_fracture_b_split.set_ylim(ymin=0, ymax=180)
    ax_particle_fracture_b_split.plot(sf_5e7.exx * 1E2, sf_5e7.fractured_binders_matrix.T,
                                      color='C1', label=sf_5e7.label)

    ax_particle_fracture_b_split.plot(sf_1e8.exx * 1E2, sf_1e8.fractured_binders_matrix.T,
                                      color='C2', label=sf_1e8.label)

    ax_particle_fracture_b_split.plot(sf_2e8.exx * 1E2, sf_2e8.fractured_binders_matrix.T,
                                      color='C3', label=sf_2e8.label)

    fig_particle_fracture_b_split.tight_layout()
    fname = fig_dir + 'particle_fracture_b_split.svg'
    plt.savefig(fname)
    ####################################################################################################################

    ##########################################BINDERN FRACTURE##########################################################
    ef_5e_2 = Simulation('ef_5e-2', no_sims, r'$\varepsilon_{f}=5$ %')
    ef_1e_1 = Simulation('ef_1e-1', no_sims, r'$\varepsilon_{f}=10$ %')
    ef_2e_1 = Simulation('ef_2e-1', no_sims, r'$\varepsilon_{f}=20$ %')

    # =FIGURE BINDER PRESSURE===========================================================================================
    fig_binder_pressure, ax_binder_pressure = plt.subplots()
    ax_binder_pressure.plot(ref.exx * 1E2,
                            ref.sxx_mean * 1E-6,
                            'k-', zorder=1)
    ax_binder_pressure.fill_between(
        ref.exx * 1E2,
        (ref.sxx_mean - ref.sxx_std) * 1E-6,
        (ref.sxx_mean + ref.sxx_std) * 1E-6,
        color='C0', label=ref.label, alpha=1, zorder=1)

    ax_binder_pressure.plot(ef_2e_1.exx * 1E2,
                            ef_2e_1.sxx_mean * 1E-6,
                            'k-', zorder=2)
    ax_binder_pressure.fill_between(
        ef_2e_1.exx * 1E2,
        (ef_2e_1.sxx_mean - ef_2e_1.sxx_std) * 1E-6,
        (ef_2e_1.sxx_mean + ef_2e_1.sxx_std) * 1E-6,
        color='C1', label=ef_2e_1.label, alpha=1, zorder=2)

    ax_binder_pressure.plot(ef_1e_1.exx * 1E2,
                            ef_1e_1.sxx_mean * 1E-6,
                            'k-', zorder=3)
    ax_binder_pressure.fill_between(
        ef_1e_1.exx * 1E2,
        (ef_1e_1.sxx_mean - ef_1e_1.sxx_std) * 1E-6,
        (ef_1e_1.sxx_mean + ef_1e_1.sxx_std) * 1E-6,
        color='C2', label=ef_1e_1.label, alpha=1, zorder=3)

    ax_binder_pressure.plot(ef_5e_2.exx * 1E2,
                            ef_5e_2.sxx_mean * 1E-6,
                            'k-', zorder=4)
    ax_binder_pressure.fill_between(
        ef_5e_2.exx * 1E2,
        (ef_5e_2.sxx_mean - ef_5e_2.sxx_std) * 1E-6,
        (ef_5e_2.sxx_mean + ef_5e_2.sxx_std) * 1E-6,
        color='C3', label=ef_5e_2.label, alpha=1, zorder=4)
    ax_binder_pressure.set_xlim(xmin=0, xmax=end_plot_strain)
    ax_binder_pressure.set_ylim(ymin=-20, ymax=120)
    ax_binder_pressure.yaxis.set_major_locator(MultipleLocator(20))
    ax_binder_pressure.xaxis.set_major_locator(MultipleLocator(1))
    ax_binder_pressure.set_ylabel('Stress [MPa]')
    ax_binder_pressure.set_xlabel('Strain [%]')


    ax_binder_pressure.legend(loc='best')
    fig_binder_pressure.tight_layout()
    fname = fig_dir + 'binder_pressure.svg'
    plt.savefig(fname)

    # =BINDER PRESSURE SPLIT============================================================================================
    fig_binder_pressure_split, ax_binder_pressure_split = plt.subplots()
    ax_binder_pressure_split.set_xlabel('Strain [%]')
    ax_binder_pressure_split.set_ylabel('Stress [MPa]')
    # ax_binder_pressure_split.set_xlim(xmin=0)
    # ax_binder_pressure_split.set_ylim(ymin=0, ymax=180)

    ax_binder_pressure_split.plot(ref.exx * 1E2, ref.in_plane_stress_matrix.T * 1E-6,
                                  color='C0', label=ref.label)

    ax_binder_pressure_split.plot(ef_5e_2.exx * 1E2, ef_5e_2.in_plane_stress_matrix.T * 1E-6,
                                  color='C1', label=ef_5e_2.label)

    ax_binder_pressure_split.plot(ef_1e_1.exx * 1E2, ef_1e_1.in_plane_stress_matrix.T * 1E-6,
                                  color='C2', label=ef_1e_1.label)

    ax_binder_pressure_split.plot(ef_2e_1.exx * 1E2, ef_2e_1.in_plane_stress_matrix.T * 1E-6,
                                  color='C3', label=ef_2e_1.label)
    ax_binder_pressure_split.legend(loc='best')
    fig_binder_pressure_split.tight_layout()
    fname = fig_dir + 'binder_pressure_split.svg'
    plt.savefig(fname)

    # =BINDER CONTACTS_P ===============================================================================================
    fig_binder_contacts_p, ax_binder_contacts_p = plt.subplots()
    ax_binder_contacts_p.plot(ref.exx * 1E2,
                              ref.particle_contacts_mean,
                              'k-', zorder=1)
    ax_binder_contacts_p.fill_between(
        ref.exx * 1E2,
        (ref.particle_contacts_mean - ref.particle_contacts_std),
        (ref.particle_contacts_mean + ref.particle_contacts_std),
        color='C0', label=ref.label, alpha=1, zorder=1)

    ax_binder_contacts_p.plot(ef_2e_1.exx * 1E2,
                              ef_2e_1.particle_contacts_mean,
                              'k-', zorder=2)
    ax_binder_contacts_p.fill_between(
        ef_2e_1.exx * 1E2,
        (ef_2e_1.particle_contacts_mean - ef_2e_1.particle_contacts_std),
        (ef_2e_1.particle_contacts_mean + ef_2e_1.particle_contacts_std),
        color='C1', label=ef_2e_1.label, alpha=1, zorder=2)

    ax_binder_contacts_p.plot(ef_1e_1.exx * 1E2,
                              ef_1e_1.particle_contacts_mean,
                              'k-', zorder=3)
    ax_binder_contacts_p.fill_between(
        ef_1e_1.exx * 1E2,
        (ef_1e_1.particle_contacts_mean - ef_1e_1.particle_contacts_std),
        (ef_1e_1.particle_contacts_mean + ef_1e_1.particle_contacts_std),
        color='C2', label=ef_1e_1.label, alpha=1, zorder=3)

    ax_binder_contacts_p.plot(ef_5e_2.exx * 1E2,
                              ef_5e_2.particle_contacts_mean,
                              'k-', zorder=4)
    ax_binder_contacts_p.fill_between(
        ef_5e_2.exx * 1E2,
        (ef_5e_2.particle_contacts_mean - ef_5e_2.particle_contacts_std),
        (ef_5e_2.particle_contacts_mean + ef_5e_2.particle_contacts_std),
        color='C3', label=ef_5e_2.label, alpha=1, zorder=4)

    ax_binder_contacts_p.set_xlim(xmin=0, xmax=end_plot_strain)
    ax_binder_contacts_p.set_ylim(ymin=2000, ymax=4000)
    ax_binder_contacts_p.yaxis.set_major_locator(MultipleLocator(250))
    ax_binder_contacts_p.xaxis.set_major_locator(MultipleLocator(1))
    ax_binder_contacts_p.set_ylabel('Particle contacts [-]')
    ax_binder_contacts_p.set_xlabel('Strain [%]')
    ax_binder_contacts_p.legend(loc='best')


    fig_binder_contacts_p.tight_layout()
    fname = fig_dir + 'binder_contacts_p.svg'
    plt.savefig(fname)

    # =BINDER CONTACTS_P SPLIT==========================================================================================
    fig_binder_contacts_p_split, ax_binder_contacts_p_split = plt.subplots()
    ax_binder_contacts_p_split.set_xlabel('Strain [%]')
    ax_binder_contacts_p_split.set_ylabel('Particle contacts [-]')
    # ax_binder_contacts_p_split.set_xlim(xmin=0)
    # ax_binder_contacts_p_split.set_ylim(ymin=0, ymax=180)
    ax_binder_contacts_p_split.plot(ref.exx * 1E2, ref.particle_contacts_matrix.T,
                                    color='C0', label=ref.label)

    ax_binder_contacts_p_split.plot(ef_5e_2.exx * 1E2, ef_5e_2.particle_contacts_matrix.T,
                                    color='C1', label=ef_5e_2.label)

    ax_binder_contacts_p_split.plot(ef_1e_1.exx * 1E2, ef_1e_1.particle_contacts_matrix.T,
                                    color='C2', label=ef_1e_1.label)

    ax_binder_contacts_p_split.plot(ef_2e_1.exx * 1E2, ef_2e_1.particle_contacts_matrix.T,
                                    color='C3', label=ef_2e_1.label)
    ax_binder_contacts_p_split.legend(loc='best')
    fig_binder_contacts_p_split.tight_layout()
    fname = fig_dir + 'binder_contacts_p_split.svg'
    plt.savefig(fname)

    # =BINDER CONTACTS_B================================================================================================
    fig_binder_contacts_b, ax_binder_contacts_b = plt.subplots()
    ax_binder_contacts_b.plot(ref.exx * 1E2,
                              ref.binder_contacts_mean,
                              'k-', zorder=1)
    ax_binder_contacts_b.fill_between(
        ref.exx * 1E2,
        (ref.binder_contacts_mean - ref.binder_contacts_std),
        (ref.binder_contacts_mean + ref.binder_contacts_std),
        color='C0', label=ref.label, alpha=1, zorder=1)

    ax_binder_contacts_b.plot(ef_2e_1.exx * 1E2,
                              ef_2e_1.binder_contacts_mean,
                              'k-', zorder=2)
    ax_binder_contacts_b.fill_between(
        ef_2e_1.exx * 1E2,
        (ef_2e_1.binder_contacts_mean - ef_2e_1.binder_contacts_std),
        (ef_2e_1.binder_contacts_mean + ef_2e_1.binder_contacts_std),
        color='C1', label=ef_2e_1.label, alpha=1, zorder=2)

    ax_binder_contacts_b.plot(ef_1e_1.exx * 1E2,
                              ef_1e_1.binder_contacts_mean,
                              'k-', zorder=3)
    ax_binder_contacts_b.fill_between(
        ef_1e_1.exx * 1E2,
        (ef_1e_1.binder_contacts_mean - ef_1e_1.binder_contacts_std),
        (ef_1e_1.binder_contacts_mean + ef_1e_1.binder_contacts_std),
        color='C2', label=ef_1e_1.label, alpha=1, zorder=3)

    ax_binder_contacts_b.plot(ef_5e_2.exx * 1E2,
                              ef_5e_2.binder_contacts_mean,
                              'k-', zorder=4)
    ax_binder_contacts_b.fill_between(
        ef_5e_2.exx * 1E2,
        (ef_5e_2.binder_contacts_mean - ef_5e_2.binder_contacts_std),
        (ef_5e_2.binder_contacts_mean + ef_5e_2.binder_contacts_std),
        color='C3', label=ef_5e_2.label, alpha=1, zorder=4)

    ax_binder_contacts_b.set_xlim(xmin=0, xmax=end_plot_strain)
    ax_binder_contacts_b.set_ylim(ymin=8750, ymax=10000)
    ax_binder_contacts_b.yaxis.set_major_locator(MultipleLocator(250))
    ax_binder_contacts_b.xaxis.set_major_locator(MultipleLocator(1))
    ax_binder_contacts_b.set_ylabel('Binder contacts [-]')
    ax_binder_contacts_b.set_xlabel('Strain [%]')


    fig_binder_contacts_b.tight_layout()
    ax_binder_contacts_b.legend(loc='best')
    fname = fig_dir + 'binder_contacts_b.svg'
    plt.savefig(fname)

    # =BINDER CONTACTS_B SPLIT==========================================================================================
    fig_binder_contacts_b_split, ax_binder_contacts_b_split = plt.subplots()
    ax_binder_contacts_b_split.set_xlabel('Strain [%]')
    ax_binder_contacts_b_split.set_ylabel('Binder contacts [-]')
    # ax_binder_contacts_b_split.set_xlim(xmin=0)
    # ax_binder_contacts_b_split.set_ylim(ymin=0, ymax=180)
    ax_binder_contacts_b_split.plot(ref.exx * 1E2, ref.binder_contacts_matrix.T,
                                    color='C0', label=ref.label)

    ax_binder_contacts_b_split.plot(ef_5e_2.exx * 1E2, ef_5e_2.binder_contacts_matrix.T,
                                    color='C1', label=ef_5e_2.label)
    ax_binder_contacts_b_split.plot(ef_1e_1.exx * 1E2, ef_1e_1.binder_contacts_matrix.T,
                                    color='C2', label=ef_1e_1.label)

    ax_binder_contacts_b_split.plot(ef_2e_1.exx * 1E2, ef_2e_1.binder_contacts_matrix.T,
                                    color='C3', label=ef_2e_1.label)
    ax_binder_contacts_b_split.legend(loc='best')
    fig_binder_contacts_b_split.tight_layout()
    fname = fig_dir + 'binder_contacts_b_split.svg'
    plt.savefig(fname)

    # =PARTICLE FRACTURE_P ===============================================================================================
    fig_binder_fracture_p, ax_binder_fracture_p = plt.subplots()
    ax_binder_fracture_p.plot(ef_2e_1.exx * 1E2,
                              ef_2e_1.fractured_particles_mean * 1E2,
                              'k-', zorder=1)
    ax_binder_fracture_p.fill_between(
        ef_2e_1.exx * 1E2,
        (ef_2e_1.fractured_particles_mean - ef_2e_1.fractured_particles_std) * 1E2,
        (ef_2e_1.fractured_particles_mean + ef_2e_1.fractured_particles_std) * 1E2,
        color='C1', label=ef_2e_1.label, alpha=1, zorder=1)

    ax_binder_fracture_p.plot(ef_1e_1.exx * 1E2,
                              ef_1e_1.fractured_particles_mean * 1E2,
                              'k-', zorder=2)
    ax_binder_fracture_p.fill_between(
        ef_1e_1.exx * 1E2,
        (ef_1e_1.fractured_particles_mean - ef_1e_1.fractured_particles_std) * 1E2,
        (ef_1e_1.fractured_particles_mean + ef_1e_1.fractured_particles_std) * 1E2,
        color='C2', label=ef_1e_1.label, alpha=1, zorder=2)

    ax_binder_fracture_p.plot(ef_5e_2.exx * 1E2,
                              ef_5e_2.fractured_particles_mean * 1E2,
                              'k-', zorder=3)
    ax_binder_fracture_p.fill_between(
        ef_5e_2.exx * 1E2,
        (ef_5e_2.fractured_particles_mean - ef_5e_2.fractured_particles_std) * 1E2,
        (ef_5e_2.fractured_particles_mean + ef_5e_2.fractured_particles_std) * 1E2,
        color='C3', label=ef_5e_2.label, alpha=1, zorder=3)

    ax_binder_fracture_p.set_xlim(xmin=0, xmax=end_plot_strain)
    ax_binder_fracture_p.set_ylim(ymin=0, ymax=100)
    ax_binder_fracture_p.yaxis.set_major_locator(MultipleLocator(10))
    ax_binder_fracture_p.xaxis.set_major_locator(MultipleLocator(1))
    ax_binder_fracture_p.set_ylabel('Fractured particles [%]')
    ax_binder_fracture_p.set_xlabel('Strain [%]')


    ax_binder_fracture_p.legend(loc='best')
    fig_binder_fracture_p.tight_layout()
    fname = fig_dir + 'binder_fracture_p.svg'
    plt.savefig(fname)

    # =BINDER FRACTURE_P SPLIT========================================================================================
    fig_binder_fracture_p_split, ax_binder_fracture_p_split = plt.subplots()
    ax_binder_fracture_p_split.set_xlabel('Strain [%]')
    ax_binder_fracture_p_split.set_ylabel('Fraction of fractured particles [-]')
    # ax_binder_fracture_p_split.set_xlim(xmin=0)
    # ax_binder_contacts_p_split.set_ylim(ymin=0, ymax=180)
    ax_binder_fracture_p_split.plot(ef_5e_2.exx * 1E2, ef_5e_2.fractured_particles_matrix.T,
                                    color='C1', label=ef_5e_2.label)
    ax_binder_fracture_p_split.plot(ef_1e_1.exx * 1E2, ef_1e_1.fractured_particles_matrix.T,
                                    color='C2', label=ef_1e_1.label)

    ax_binder_fracture_p_split.plot(ef_2e_1.exx * 1E2, ef_2e_1.fractured_particles_matrix.T,
                                    color='C3', label=ef_2e_1.label)

    ax_binder_fracture_p_split.legend(loc='best')
    fig_binder_fracture_p_split.tight_layout()
    fname = fig_dir + 'binder_fracture_p_split.svg'
    plt.savefig(fname)

    # =BINDER FRACTURE_B================================================================================================
    fig_binder_fracture_b, ax_binder_fracture_b = plt.subplots()
    ax_binder_fracture_b.plot(ef_2e_1.exx * 1E2,
                              ef_2e_1.fractured_binders_mean * 1E2,
                              'k-', zorder=1)
    ax_binder_fracture_b.fill_between(
        ef_2e_1.exx * 1E2,
        (ef_2e_1.fractured_binders_mean - ef_2e_1.fractured_binders_std) * 1E2,
        (ef_2e_1.fractured_binders_mean + ef_2e_1.fractured_binders_std) * 1E2,
        color='C1', label=ef_2e_1.label, alpha=1, zorder=1)

    ax_binder_fracture_b.plot(ef_1e_1.exx * 1E2,
                              ef_1e_1.fractured_binders_mean * 1E2,
                              'k-', zorder=2)
    ax_binder_fracture_b.fill_between(
        ef_1e_1.exx * 1E2,
        (ef_1e_1.fractured_binders_mean - ef_1e_1.fractured_binders_std) * 1E2,
        (ef_1e_1.fractured_binders_mean + ef_1e_1.fractured_binders_std) * 1E2,
        color='C2', label=ef_1e_1.label, alpha=1, zorder=2)

    ax_binder_fracture_b.plot(ef_5e_2.exx * 1E2,
                              ef_5e_2.fractured_binders_mean * 1E2,
                              'k-', zorder=3)
    ax_binder_fracture_b.fill_between(
        ef_5e_2.exx * 1E2,
        (ef_5e_2.fractured_binders_mean - ef_5e_2.fractured_binders_std) * 1E2,
        (ef_5e_2.fractured_binders_mean + ef_5e_2.fractured_binders_std) * 1E2,
        color='C3', label=ef_5e_2.label, alpha=1, zorder=3)

    ax_binder_fracture_b.set_xlim(xmin=0, xmax=end_plot_strain)
    ax_binder_fracture_b.set_ylim(ymin=0, ymax=25)
    ax_binder_fracture_b.yaxis.set_major_locator(MultipleLocator(5))
    ax_binder_fracture_b.xaxis.set_major_locator(MultipleLocator(1))
    ax_binder_fracture_b.set_ylabel('Fractured binder contacts [%]')
    ax_binder_fracture_b.set_xlabel('Strain [%]')

    fig_binder_fracture_b.tight_layout()
    ax_binder_fracture_b.legend(loc='best')
    fname = fig_dir + 'binder_fracture_b.svg'
    plt.savefig(fname)

    # =PARTICLE FRACTURE_B SPLIT========================================================================================
    fig_binder_fracture_b_split, ax_binder_fracture_b_split = plt.subplots()
    ax_binder_fracture_b_split.set_xlabel('Strain [%]')
    ax_binder_fracture_b_split.set_ylabel('Fraction of fractured binder contacts [-]')
    # ax_binder_fracture_b_split.set_xlim(xmin=0)
    # ax_binder_fracture_b_split.set_ylim(ymin=0, ymax=180)
    ax_binder_fracture_b_split.plot(ef_5e_2.exx * 1E2, ef_5e_2.fractured_binders_matrix.T,
                                    color='C1', label=ef_5e_2.label)

    ax_binder_fracture_b_split.plot(ef_1e_1.exx * 1E2, ef_1e_1.fractured_binders_matrix.T,
                                    color='C2', label=ef_1e_1.label)

    ax_binder_fracture_b_split.plot(ef_2e_1.exx * 1E2, ef_2e_1.fractured_binders_matrix.T,
                                    color='C3', label=ef_2e_1.label)

    fig_binder_fracture_b_split.tight_layout()
    fname = fig_dir + 'binder_fracture_b_split.svg'
    plt.savefig(fname)
    ####################################################################################################################

    #######################################BINDER AND PARTICLE FRACTURE#################################################

    sf_2e8_ef_1e_1 = Simulation('sf_2e8_ef_1e-1', no_sims, r'$\varepsilon=0$ %')

    # =FIGURE BINDER PARTICLE PRESSURE==================================================================================
    fig_particle_binder_pressure, ax_particle_binder_pressure = plt.subplots()
    ax_particle_binder_pressure.plot(sf_2e8_ef_1e_1.exx * 1E2,
                                     sf_2e8_ef_1e_1.sxx_mean* 1E-6,
                                     'k-', zorder=1)
    ax_particle_binder_pressure.fill_between(
        sf_2e8_ef_1e_1.exx * 1E2,
        (sf_2e8_ef_1e_1.sxx_mean - sf_2e8_ef_1e_1.sxx_std) * 1E-6,
        (sf_2e8_ef_1e_1.sxx_mean + sf_2e8_ef_1e_1.sxx_std) * 1E-6,
        color='C2', label=r'$\sigma_{xx}$', alpha=1, zorder=1)

    ax_particle_binder_pressure.set_xlim(xmin=0, xmax=end_plot_strain)
    ax_particle_binder_pressure.set_ylim(ymin=-10, ymax=40)
    ax_particle_binder_pressure.yaxis.set_major_locator(MultipleLocator(10))
    ax_particle_binder_pressure.xaxis.set_major_locator(MultipleLocator(1))
    ax_particle_binder_pressure.set_ylabel('Stress [MPa]')
    ax_particle_binder_pressure.set_xlabel('Strain [%]')

    fig_particle_binder_pressure.tight_layout()
    # ax_particle_binder_pressure.legend(loc='best')
    fname = fig_dir + 'particle_binder_pressure.svg'
    plt.savefig(fname)

    # =BINDER PARTICLE CONTACTS_P ======================================================================================
    fig_particle_binder_contacts_p, ax_particle_binder_contacts_p = plt.subplots()
    ax_particle_binder_contacts_p.plot(sf_2e8_ef_1e_1.exx * 1E2,
                                       sf_2e8_ef_1e_1.particle_contacts_mean,
                                       'k-', zorder=2)
    ax_particle_binder_contacts_p.fill_between(
        sf_2e8_ef_1e_1.exx * 1E2,
        (sf_2e8_ef_1e_1.particle_contacts_mean - sf_2e8_ef_1e_1.particle_contacts_std),
        (sf_2e8_ef_1e_1.particle_contacts_mean + sf_2e8_ef_1e_1.particle_contacts_std),
        color='C2', label=sf_2e8_ef_1e_1.label, alpha=1, zorder=2)

    ax_particle_binder_contacts_p.set_xlim(xmin=0,xmax=end_plot_strain)
    ax_particle_binder_contacts_p.set_ylim(ymin=2400, ymax=3400)
    ax_particle_binder_contacts_p.xaxis.set_major_locator(MultipleLocator(1))
    ax_particle_binder_contacts_p.yaxis.set_major_locator(MultipleLocator(200))
    ax_particle_binder_contacts_p.set_ylabel('Particle contacts [-]')
    ax_particle_binder_contacts_p.set_xlabel('Strain [%]')

    # ax_particle_binder_contacts_p.legend(loc='best')
    fig_particle_binder_contacts_p.tight_layout()
    fname = fig_dir + 'particle_binder_contacts_p.svg'
    plt.savefig(fname)

    # =BINDER PARTICLE CONTACTS_P SPLIT=================================================================================
    fig_particle_binder_contacts_p_split, ax_particle_binder_contacts_p_split = plt.subplots()
    ax_particle_binder_contacts_p_split.set_xlabel('Strain [%]')
    ax_particle_binder_contacts_p_split.set_ylabel('Particle contacts [-]')

    ax_particle_binder_contacts_p_split.plot(sf_2e8_ef_1e_1.exx * 1E2,
                                             sf_2e8_ef_1e_1.particle_contacts_matrix.T,
                                             color='C2', label=sf_2e8_ef_1e_1.label)

    ax_particle_binder_contacts_p_split.legend(loc='best')
    fig_particle_binder_contacts_p_split.tight_layout()
    fname = fig_dir + 'particle_binder_contacts_p_split.svg'
    plt.savefig(fname)

    # =BINDER PARTICLE CONTACTS_B=======================================================================================
    fig_particle_binder_contacts_b, ax_particle_binder_contacts_b = plt.subplots()
    ax_particle_binder_contacts_b.plot(sf_2e8_ef_1e_1.exx * 1E2,
                                       sf_2e8_ef_1e_1.binder_contacts_mean,
                                       'k-', zorder=2)
    ax_particle_binder_contacts_b.fill_between(
        sf_2e8_ef_1e_1.exx * 1E2,
        (sf_2e8_ef_1e_1.binder_contacts_mean - sf_2e8_ef_1e_1.binder_contacts_std),
        (sf_2e8_ef_1e_1.binder_contacts_mean + sf_2e8_ef_1e_1.binder_contacts_std),
        color='C2', label=sf_2e8_ef_1e_1.label, alpha=1, zorder=2)

    ax_particle_binder_contacts_b.set_xlim(xmin=0,xmax=end_plot_strain)
    ax_particle_binder_contacts_b.set_ylim(ymin=9200, ymax=10000)
    ax_particle_binder_contacts_b.xaxis.set_major_locator(MultipleLocator(1))
    ax_particle_binder_contacts_b.yaxis.set_major_locator(MultipleLocator(100))
    ax_particle_binder_contacts_b.set_ylabel('Binder contacts [-]')
    ax_particle_binder_contacts_b.set_xlabel('Strain [%]')


    fig_particle_binder_contacts_b.tight_layout()
    # ax_particle_binder_contacts_b.legend(loc='best')
    fname = fig_dir + 'particle_binder_contacts_b.svg'
    plt.savefig(fname)

    # =BINDER PARTICLE CONTACTS_B SPLIT=================================================================================
    fig_particle_binder_contacts_b_split, ax_particle_binder_contacts_b_split = plt.subplots()
    ax_particle_binder_contacts_b_split.set_xlabel('Strain [%]')
    ax_particle_binder_contacts_b_split.set_ylabel('Binder contacts [-]')
    # ax_particle_binder_contacts_b_split.set_xlim(xmin=0)
    # ax_particle_binder_contacts_b_split.set_ylim(ymin=0, ymax=180)

    ax_particle_binder_contacts_b_split.plot(sf_2e8_ef_1e_1.exx * 1E2,
                                             sf_2e8_ef_1e_1.binder_contacts_matrix.T,
                                             color='C2', label=sf_2e8_ef_1e_1.label)


    ax_particle_binder_contacts_b_split.legend(loc='best')
    fig_particle_binder_contacts_b_split.tight_layout()
    fname = fig_dir + 'particle_binder_contacts_b_split.svg'
    plt.savefig(fname)

    # =BINDER PARTICLE FRACTURE_P ======================================================================================
    fig_particle_binder_fracture_p, ax_particle_binder_fracture_p = plt.subplots()
    ax_particle_binder_fracture_p.plot(sf_2e8_ef_1e_1.exx * 1E2,
                                       sf_2e8_ef_1e_1.fractured_particles_mean,
                                       'k-', zorder=2)
    ax_particle_binder_fracture_p.fill_between(
        sf_2e8_ef_1e_1.exx * 1E2,
        (sf_2e8_ef_1e_1.fractured_particles_mean - sf_2e8_ef_1e_1.fractured_particles_std),
        (sf_2e8_ef_1e_1.fractured_particles_mean + sf_2e8_ef_1e_1.fractured_particles_std),
        color='C2', label=sf_2e8_ef_1e_1.label, alpha=1, zorder=2)

    ax_particle_binder_fracture_p.set_xlim(xmin=0,xmax=end_plot_strain)
    ax_particle_binder_fracture_p.xaxis.set_major_locator(MultipleLocator(1))
    ax_particle_binder_fracture_p.set_ylim(ymin=0, ymax=1.0)
    ax_particle_binder_fracture_p.yaxis.set_major_locator(MultipleLocator(0.1))
    ax_particle_binder_fracture_p.set_ylabel('Fractured particles [%]')
    ax_particle_binder_fracture_p.set_xlabel('Strain [%]')

    # ax_particle_binder_fracture_p.legend(loc='best')
    fig_particle_binder_fracture_p.tight_layout()
    fname = fig_dir + 'particle_binder_fracture_p.svg'
    plt.savefig(fname)

    # =BINDER PARTICLE FRACTURE_P SPLIT=================================================================================
    fig_particle_binder_fracture_p_split, ax_particle_binder_fracture_p_split = plt.subplots()
    ax_particle_binder_fracture_p_split.set_xlabel('Strain [%]')
    ax_particle_binder_fracture_p_split.set_ylabel('Fraction of fractured particles [-]')
    ax_particle_binder_fracture_p_split.plot(sf_2e8_ef_1e_1.exx * 1E2,
                                             sf_2e8_ef_1e_1.fractured_particles_matrix.T,
                                             color='C2', label=sf_2e8_ef_1e_1.label)

    ax_particle_binder_fracture_p_split.legend(loc='best')
    fig_particle_binder_fracture_p_split.tight_layout()
    fname = fig_dir + 'particle_binder_fracture_p_split.svg'
    plt.savefig(fname)

    # =BINDER PARTICLE FRACTURE_B=======================================================================================
    fig_particle_binder_fracture_b, ax_particle_binder_fracture_b = plt.subplots()
    ax_particle_binder_fracture_b.plot(sf_2e8_ef_1e_1.exx * 1E2,
                                       sf_2e8_ef_1e_1.fractured_binders_mean,
                                       'k-', zorder=2)
    ax_particle_binder_fracture_b.fill_between(
        sf_2e8_ef_1e_1.exx * 1E2,
        (sf_2e8_ef_1e_1.fractured_binders_mean - sf_2e8_ef_1e_1.fractured_binders_std),
        (sf_2e8_ef_1e_1.fractured_binders_mean + sf_2e8_ef_1e_1.fractured_binders_std),
        color='C2', label=sf_2e8_ef_1e_1.label, alpha=1, zorder=2)

    ax_particle_binder_fracture_b.set_xlim(xmin=0,xmax=end_plot_strain)
    ax_particle_binder_fracture_b.xaxis.set_major_locator(MultipleLocator(1))
    ax_particle_binder_fracture_b.set_ylim(ymin=0, ymax=0.10)
    ax_particle_binder_fracture_b.yaxis.set_major_locator(MultipleLocator(0.01))
    ax_particle_binder_fracture_b.set_ylabel('Fractured binder contacts [%]')
    ax_particle_binder_fracture_b.set_xlabel('Strain [%]')

    # ax_particle_binder_fracture_b.legend(loc='best')
    fig_particle_binder_fracture_b.tight_layout()
    fname = fig_dir + 'particle_binder_fracture_b.svg'
    plt.savefig(fname)

    # =BINDER PARTICLE FRACTURE_B SPLIT=================================================================================
    fig_particle_binder_fracture_b_split, ax_particle_binder_fracture_b_split = plt.subplots()
    ax_particle_binder_fracture_b_split.set_xlabel('Strain [%]')
    ax_particle_binder_fracture_b_split.set_ylabel('Fraction of fractured binder contacts [-]')
    ax_particle_binder_fracture_b_split.plot(sf_2e8_ef_1e_1.exx * 1E2,
                                             sf_2e8_ef_1e_1.fractured_binders_matrix.T,
                                             color='C2', label=sf_2e8_ef_1e_1.label)

    fig_particle_binder_fracture_b_split.tight_layout()
    fname = fig_dir + 'particle_binder_fracture_b_split.svg'
    plt.savefig(fname)

    ########################################################################################################################
    plt.show()
