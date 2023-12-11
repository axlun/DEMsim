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


def stiffness_spread_processing(sim_dir, sim_type, no_sims):
    simulation_directory = {}

    for i in range(1, no_sims + 1):
        simulation_directory[i] = MechanicalLoading(sim_dir + str(i) + '/electrode_mechanical_loading_' + str(sim_type))

    stiffness_points_matrix = np.zeros((no_sims, simulation_directory[1].strain_points_total.size))
    strain_points = simulation_directory[1].strain_points_total
    for i in range(no_sims):
        stiffness_points_matrix[i, :] = simulation_directory[i + 1].stiffness_values_total

    stiffness_mean = np.mean(stiffness_points_matrix, axis=0)
    stiffness_std = np.std(stiffness_points_matrix, axis=0)
    return strain_points, stiffness_mean, stiffness_std


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
        self.sxx, self.syy, self.szz, self.tau_xy, self.tau_xz, self.tau_yx, self.tau_yz, self.tau_zx, self.tau_zy, \
        self.ke = calendering_plot_processing(sim_dir)
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


class MechanicalLoading:
    def __init__(self, sim_dir):
        self.sim_dir = sim_dir

        self.time_tension, self.linear_strain_tension, self.sxx_tension, self.syy_tension, self.szz_tension, \
            self.tau_xy_tension, self.tau_xz_tension, self.tau_yx_tension, self.tau_yz_tension, self.tau_zx_tension, \
            self.tau_zy_tension, self.time_compression, self.linear_strain_compression, self.sxx_compression, \
            self.syy_compression, self.szz_compression, self.tau_xy_compression, self.tau_xz, self.tau_yx, \
            self.tau_yz_compression, self.tau_zx_compression, self.tau_zy_compression = mech_plot_prop(sim_dir)
        self.strain_points_total, self.stiffness_values_total = stiffness_func(self.sxx_tension,
                                                                               self.linear_strain_tension,
                                                                               self.sxx_compression,
                                                                               self.linear_strain_compression)


class Simulation:
    def __init__(self, sim_dir, sim_type, no_sims, label='template'):
        self.sim_dir = sim_dir
        self.sim_type = sim_type
        self.no_sims = no_sims
        self.strain_points, self.stiffness_mean, self.stiffness_std = stiffness_spread_processing(sim_dir, sim_type,
                                                                                                  no_sims)
        self.loading_surface_position, self.loading_surface_pressure_mean, self.loading_surface_pressure_std, \
            self.unloading_surface_position, self.unloading_surface_pressure_mean, self.unloading_surface_pressure_std \
            = calendering_spread_processing(sim_dir, sim_type, no_sims)
        self.label = label


if __name__ == '__main__':

    simulation_directory_SN_101 = '/scratch/users/axlun/DEMsim/results/article_2/final_runs_2/SN_101/'
    SN_101 = Simulation(simulation_directory_SN_101, 'hertz', 4, 'SN_101')

    simulation_directory_SN_101P = '/scratch/users/axlun/DEMsim/results/article_2/final_runs_2/SN_101P/'
    SN_101P = Simulation(simulation_directory_SN_101P, 'hertz', 4, 'SN_101P')

    simulation_directory_SN_111 = '/scratch/users/axlun/DEMsim/results/article_2/final_runs_2/SN_111/'
    SN_111 = Simulation(simulation_directory_SN_111, 'hertz', 4, 'SN_111')

    simulation_directory_SN_200 = '/scratch/users/axlun/DEMsim/results/article_2/final_runs_2/SN_200/'
    SN_200 = Simulation(simulation_directory_SN_200, 'el_pl_binder_el_pl_particle', 4, 'SN_200')

    simulation_directory_SN_201 = '/scratch/users/axlun/DEMsim/results/article_2/final_runs_2/SN_201/'
    SN_201 = Simulation(simulation_directory_SN_201, 'el_pl_binder_el_pl_particle', 4, 'SN_201')

    # simulation_directory_SN_202 = '/scratch/users/axlun/DEMsim/results/article_2/final_runs_2/SN_202/'
    # SN_202 = Simulation(simulation_directory_SN_202, 'el_pl_binder_el_pl_particle', 4, 'SN_202')

    simulation_directory_SN_301 = '/scratch/users/axlun/DEMsim/results/article_2/final_runs_2/SN_301/'
    SN_301 = Simulation(simulation_directory_SN_301, 'el_pl_binder_el_pl_particle', 4, 'SN_301')

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

    # ==EXPERIMENTAL DATA===============================================================================================
    exp_strain_points_compression = [-1.00, -1.10, -1.23, -1.41, -1.65]
    """
    #Wrong values read from from experiment result table
    Modulus_eps_dot_01_compression = [1.20, 1.04, 0.96, 1.19, 1.46]
    Modulus_eps_dot_05_compression = [1.43, 0.86, 0.97, 1.39, 1.79]
    Modulus_eps_dot_10_compression = [1.51, 1.79, 1.99, 2.09, 2.74]
    Modulus_eps_dot_100_compression = [1.55, 1.72, 1.70, 1.81, 2.28]
    Modulus_eps_dot_300_compression = [1.99, 1.52, 2.21, 2.80, 2.33]
    """
    Modulus_eps_dot_01_compression = [1.20, 1.43, 1.51, 1.55, 1.99]
    Modulus_eps_dot_05_compression = [1.04, 0.86, 1.79, 1.72, 1.52]
    Modulus_eps_dot_10_compression = [0.96, 0.97, 1.99, 1.70, 2.21]
    Modulus_eps_dot_100_compression = [1.19, 1.39, 2.09, 1.81, 2.80]
    Modulus_eps_dot_300_compression = [1.46, 1.79, 2.74, 2.28, 2.33]

    exp_strain_points_tension = [1.00, 1.10, 1.23, 1.41, 1.65]
    """
    #Wrong values read from from experiment result table
    Modulus_eps_dot_01_tension = [0.95, 0.50, 0.89, 0.61, 0.57]
    Modulus_eps_dot_05_tension = [0.78, 0.39, 1.17, 0.55, 0.46]
    Modulus_eps_dot_10_tension = [0.84, 0.59, 1.09, 0.59, 0.52]
    Modulus_eps_dot_100_tension = [0.76, 0.55, 0.91, 0.63, 0.67]
    Modulus_eps_dot_300_tension = [1.06, 0.55, 0.90, 0.95, 0.85]
    """
    Modulus_eps_dot_01_tension = [0.95, 0.78, 0.84, 0.76, 1.06]
    Modulus_eps_dot_05_tension = [0.50, 0.39, 0.59, 0.55, 0.55]
    Modulus_eps_dot_10_tension = [0.89, 1.17, 1.09, 0.91, 0.90]
    Modulus_eps_dot_100_tension = [0.61, 0.55, 0.59, 0.63, 0.95]
    Modulus_eps_dot_300_tension = [0.57, 0.46, 0.52, 0.67, 0.85]

    exp_strain_points = exp_strain_points_compression[::-1] + exp_strain_points_tension
    Modulus_eps_dot_01 = Modulus_eps_dot_01_compression[::-1] + Modulus_eps_dot_01_tension
    Modulus_eps_dot_05 = Modulus_eps_dot_05_compression[::-1] + Modulus_eps_dot_05_tension
    Modulus_eps_dot_10 = Modulus_eps_dot_10_compression[::-1] + Modulus_eps_dot_10_tension
    Modulus_eps_dot_100 = Modulus_eps_dot_100_compression[::-1] + Modulus_eps_dot_100_tension
    Modulus_eps_dot_300 = Modulus_eps_dot_300_compression[::-1] + Modulus_eps_dot_300_tension
    # =FIGURES==========================================================================================================

    # =FIGURE PARTICLE SIZE=============================================================================================
    fig_PSD, ax_PSD = plt.subplots()
    ax_PSD.set_ylabel('Unloading stiffness [GPa]')
    ax_PSD.set_xlabel('Strain [%]')
    # --SIMULATIONS----------------------------------------------------------------------------------------------------
    lns_SN_101_mean = ax_PSD.plot(SN_101.strain_points * 1E2, SN_101.stiffness_mean * 1E-9,
                                  linestyle='dashed', marker='x', markersize=12, markeredgewidth=3, linewidth=3,
                                  color='C0', label=SN_101.label)
    lns_SN_101_std = ax_PSD.errorbar(SN_101.strain_points * 1E2, SN_101.stiffness_mean * 1E-9,
                                     yerr=SN_101.stiffness_std * 1E-9, elinewidth=2, capthick=None, capsize=12,
                                     lw=0, marker='', markersize=0, markeredgewidth=2, color='C0',
                                     label='Standard deviation')
    lns_SN_101P_mean = ax_PSD.plot(SN_101P.strain_points * 1E2, SN_101P.stiffness_mean * 1E-9,
                                   linestyle='dashed', marker='x', markersize=12, markeredgewidth=3, linewidth=3,
                                   color='C1', label=SN_101P.label)
    lns_SN_101P_std = ax_PSD.errorbar(SN_101P.strain_points * 1E2, SN_101P.stiffness_mean * 1E-9,
                                      yerr=SN_101P.stiffness_std * 1E-9, elinewidth=2, capthick=None, capsize=12,
                                      lw=0, marker='', markersize=0, markeredgewidth=2, color='C1',
                                      label='Standard deviation')

    # --EXPERIMENTS-----------------------------------------------------------------------------------------------------
    lns_stiff_exp_eps_dot_01_El_Pl = ax_PSD.plot(exp_strain_points, Modulus_eps_dot_01,
                                                 marker='o', markersize=10, markeredgewidth=3,
                                                 linewidth=0, color='C9',
                                                 label=r'$\dot{\Delta} = 0.1 mm/min$')
    lns_stiff_exp_eps_dot_05_El_Pl = ax_PSD.plot(exp_strain_points, Modulus_eps_dot_05,
                                                 marker='v', markersize=10, markeredgewidth=3,
                                                 linewidth=0, color='C8',
                                                 label='$\dot{\Delta} = 0.5 mm/min$')
    lns_stiff_exp_eps_dot_10_El_Pl = ax_PSD.plot(exp_strain_points, Modulus_eps_dot_10,
                                                 marker='s', markersize=10, markeredgewidth=3,
                                                 linewidth=0, color='C7',
                                                 label='$\dot{\Delta} = 1.0 mm/min$')
    lns_stiff_exp_eps_dot_100_El_Pl = ax_PSD.plot(exp_strain_points, Modulus_eps_dot_100,
                                                  marker='P', markersize=10, markeredgewidth=3,
                                                  linewidth=0, color='C6',
                                                  label='$\dot{\Delta} = 10 mm/min$')
    lns_stiff_exp_eps_dot_300_El_Pl = ax_PSD.plot(exp_strain_points, Modulus_eps_dot_300,
                                                  marker='D', markersize=10, markeredgewidth=3,
                                                  linewidth=0, color='C5',
                                                  label='$\dot{\Delta} = 30 mm/min$')
    # ===========Legends=================================
    handles_exp = lns_stiff_exp_eps_dot_01_El_Pl + lns_stiff_exp_eps_dot_05_El_Pl + lns_stiff_exp_eps_dot_10_El_Pl + \
                  lns_stiff_exp_eps_dot_100_El_Pl + lns_stiff_exp_eps_dot_300_El_Pl
    handles_exp.reverse()
    labels_exp = [l.get_label() for l in handles_exp]
    first_legend_exp = plt.legend(handles_exp, labels_exp, loc='upper right', title='Experimental results')

    ax_PSD.add_artist(first_legend_exp)

    handles_sim = lns_SN_101_mean + lns_SN_101P_mean
    labels_sim = [l.get_label() for l in handles_sim]
    plt.legend(handles_sim, labels_sim, loc='upper center', title='Simulations')

    ax_PSD.set_ylim(ymin=0)
    fname = fig_dir + 'particle_size'
    plt.savefig(fname)

    # =FIGURE CONTACT MODEL=============================================================================================
    fig_stiffness, ax_stiffness = plt.subplots()
    ax_stiffness.set_ylabel('Unloading stiffness [GPa]')
    ax_stiffness.set_xlabel('Strain [%]')
    # ax_stiff_El_Pl.set_title('Unloading stiffness of electrode layer')
    # --SIMULATIONS----------------------------------------------------------------------------------------------------
    lns_SN_101_mean = ax_stiffness.plot(SN_101.strain_points * 1E2, SN_101.stiffness_mean * 1E-9,
                                        linestyle='dashed', marker='x', markersize=12, markeredgewidth=3, linewidth=3,
                                        color='C0', label=SN_101.label)
    lns_SN_101_std = ax_stiffness.errorbar(SN_101.strain_points * 1E2, SN_101.stiffness_mean * 1E-9,
                                           yerr=SN_101.stiffness_std * 1E-9, elinewidth=2, capthick=None, capsize=12,
                                           lw=0, marker='', markersize=0, markeredgewidth=2, color='C0',
                                           label='Standard deviation')

    # lns_SN_111_mean= ax_stiffness.plot(SN_111.strain_points* 1E2, SN_111.stiffness_mean* 1E-9,
    #                                   linestyle='dashed', marker='x', markersize=12, markeredgewidth=3, linewidth=3,
    #                                   color='C1', label=SN_111.label)
    # lns_SN_111_std = ax_stiffness.errorbar(SN_111.strain_points * 1E2, SN_111.stiffness_mean*1E-9,
    #                                       yerr=SN_111.stiffness_std*1E-9, elinewidth=2, capthick=None, capsize=12,
    #                                       lw=0, marker='', markersize=0, markeredgewidth=2, color='C1',
    #                                       label='Standard deviation')

    lns_SN_201_mean = ax_stiffness.plot(SN_201.strain_points * 1E2, SN_201.stiffness_mean * 1E-9,
                                        linestyle='dashed', marker='x', markersize=12, markeredgewidth=3, linewidth=3,
                                        color='C2', label=SN_201.label)
    lns_SN_201_std = ax_stiffness.errorbar(SN_201.strain_points * 1E2, SN_201.stiffness_mean * 1E-9,
                                           yerr=SN_201.stiffness_std * 1E-9, elinewidth=2, capthick=None, capsize=12,
                                           lw=0, marker='', markersize=0, markeredgewidth=2, color='C2',
                                           label='Standard deviation')
    lns_SN_301_mean = ax_stiffness.plot(SN_301.strain_points * 1E2, SN_301.stiffness_mean * 1E-9,
                                        linestyle='dashed', marker='x', markersize=12, markeredgewidth=3, linewidth=3,
                                        color='C3', label=SN_301.label)
    lns_SN_301_std = ax_stiffness.errorbar(SN_301.strain_points * 1E2, SN_301.stiffness_mean * 1E-9,
                                           yerr=SN_301.stiffness_std * 1E-9, elinewidth=2, capthick=None, capsize=12,
                                           lw=0, marker='', markersize=0, markeredgewidth=2, color='C3',
                                           label='Standard deviation')
    # --EXPERIMENTS-----------------------------------------------------------------------------------------------------
    lns_stiff_exp_eps_dot_01_El_Pl = ax_stiffness.plot(exp_strain_points, Modulus_eps_dot_01,
                                                       marker='o', markersize=10, markeredgewidth=3,
                                                       linewidth=0, color='C9',
                                                       label=r'$\dot{\Delta} = 0.1 mm/min$')
    lns_stiff_exp_eps_dot_05_El_Pl = ax_stiffness.plot(exp_strain_points, Modulus_eps_dot_05,
                                                       marker='v', markersize=10, markeredgewidth=3,
                                                       linewidth=0, color='C8',
                                                       label='$\dot{\Delta} = 0.5 mm/min$')
    lns_stiff_exp_eps_dot_10_El_Pl = ax_stiffness.plot(exp_strain_points, Modulus_eps_dot_10,
                                                       marker='s', markersize=10, markeredgewidth=3,
                                                       linewidth=0, color='C7',
                                                       label='$\dot{\Delta} = 1.0 mm/min$')
    lns_stiff_exp_eps_dot_100_El_Pl = ax_stiffness.plot(exp_strain_points, Modulus_eps_dot_100,
                                                        marker='P', markersize=10, markeredgewidth=3,
                                                        linewidth=0, color='C6',
                                                        label='$\dot{\Delta} = 10 mm/min$')
    lns_stiff_exp_eps_dot_300_El_Pl = ax_stiffness.plot(exp_strain_points, Modulus_eps_dot_300,
                                                        marker='D', markersize=10, markeredgewidth=3,
                                                        linewidth=0, color='C5',
                                                        label='$\dot{\Delta} = 30 mm/min$')
    # ===========Legends=================================
    handles_exp = lns_stiff_exp_eps_dot_01_El_Pl + lns_stiff_exp_eps_dot_05_El_Pl + lns_stiff_exp_eps_dot_10_El_Pl + \
                  lns_stiff_exp_eps_dot_100_El_Pl + lns_stiff_exp_eps_dot_300_El_Pl
    handles_exp.reverse()
    labels_exp = [l.get_label() for l in handles_exp]
    first_legend_exp = plt.legend(handles_exp, labels_exp, loc='upper right', title='Experimental results')

    ax_stiffness.add_artist(first_legend_exp)

    handles_sim = lns_SN_101_mean + lns_SN_201_mean + lns_SN_301_mean
    labels_sim = [l.get_label() for l in handles_sim]
    plt.legend(handles_sim, labels_sim, loc='upper center', title='Simulations')

    ax_stiffness.set_ylim(ymin=0)
    fname = fig_dir + 'contact_model'
    plt.savefig(fname)

    # =FIGURE CALEDNERING DEGREE========================================================================================
    fig_cal_degree, ax_cal_degree = plt.subplots()
    ax_cal_degree.set_ylabel('Unloading stiffness [GPa]')
    ax_cal_degree.set_xlabel('Strain [%]')
    # --SIMULATIONS----------------------------------------------------------------------------------------------------
    lns_SN_200_mean = ax_cal_degree.plot(SN_200.strain_points * 1E2, SN_200.stiffness_mean * 1E-9,
                                         linestyle='dashed', marker='x', markersize=12, markeredgewidth=3, linewidth=3,
                                         color='C0', label=SN_200.label)
    lns_SN_200_std = ax_cal_degree.errorbar(SN_200.strain_points * 1E2, SN_200.stiffness_mean * 1E-9,
                                            yerr=SN_200.stiffness_std * 1E-9, elinewidth=2, capthick=None, capsize=12,
                                            lw=0, marker='', markersize=0, markeredgewidth=2, color='C0',
                                            label='Standard deviation')

    lns_SN_201_mean = ax_cal_degree.plot(SN_201.strain_points * 1E2, SN_201.stiffness_mean * 1E-9,
                                         linestyle='dashed', marker='x', markersize=12, markeredgewidth=3, linewidth=3,
                                         color='C2', label=SN_201.label)
    lns_SN_201_std = ax_cal_degree.errorbar(SN_201.strain_points * 1E2, SN_201.stiffness_mean * 1E-9,
                                            yerr=SN_201.stiffness_std * 1E-9, elinewidth=2, capthick=None, capsize=12,
                                            lw=0, marker='', markersize=0, markeredgewidth=2, color='C2',
                                            label='Standard deviation')
    """
    lns_SN_202_mean = ax_cal_degree.plot(SN_202.strain_points * 1E2, SN_202.stiffness_mean * 1E-9,
                                         linestyle='dashed', marker='x', markersize=12, markeredgewidth=3, linewidth=3,
                                         color='C3', label=SN_202.label)
    lns_SN_202_std = ax_cal_degree.errorbar(SN_202.strain_points * 1E2, SN_202.stiffness_mean * 1E-9,
                                            yerr=SN_202.stiffness_std * 1E-9, elinewidth=2, capthick=None, capsize=12,
                                            lw=0, marker='', markersize=0, markeredgewidth=2, color='C3',
                                            label='Standard deviation')
    """
    # --EXPERIMENTS-----------------------------------------------------------------------------------------------------
    lns_stiff_exp_eps_dot_01_El_Pl = ax_cal_degree.plot(exp_strain_points, Modulus_eps_dot_01,
                                                        marker='o', markersize=10, markeredgewidth=3,
                                                        linewidth=0, color='C9',
                                                        label=r'$\dot{\Delta} = 0.1 mm/min$')
    lns_stiff_exp_eps_dot_05_El_Pl = ax_cal_degree.plot(exp_strain_points, Modulus_eps_dot_05,
                                                        marker='v', markersize=10, markeredgewidth=3,
                                                        linewidth=0, color='C8',
                                                        label='$\dot{\Delta} = 0.5 mm/min$')
    lns_stiff_exp_eps_dot_10_El_Pl = ax_cal_degree.plot(exp_strain_points, Modulus_eps_dot_10,
                                                        marker='s', markersize=10, markeredgewidth=3,
                                                        linewidth=0, color='C7',
                                                        label='$\dot{\Delta} = 1.0 mm/min$')
    lns_stiff_exp_eps_dot_100_El_Pl = ax_cal_degree.plot(exp_strain_points, Modulus_eps_dot_100,
                                                         marker='P', markersize=10, markeredgewidth=3,
                                                         linewidth=0, color='C6',
                                                         label='$\dot{\Delta} = 10 mm/min$')
    lns_stiff_exp_eps_dot_300_El_Pl = ax_cal_degree.plot(exp_strain_points, Modulus_eps_dot_300,
                                                         marker='D', markersize=10, markeredgewidth=3,
                                                         linewidth=0, color='C5',
                                                         label='$\dot{\Delta} = 30 mm/min$')
    # ===========Legends=================================
    handles_exp = lns_stiff_exp_eps_dot_01_El_Pl + lns_stiff_exp_eps_dot_05_El_Pl + lns_stiff_exp_eps_dot_10_El_Pl + \
                  lns_stiff_exp_eps_dot_100_El_Pl + lns_stiff_exp_eps_dot_300_El_Pl
    handles_exp.reverse()
    labels_exp = [l.get_label() for l in handles_exp]
    first_legend_exp = plt.legend(handles_exp, labels_exp, loc='upper right', title='Experimental results')

    ax_cal_degree.add_artist(first_legend_exp)

    handles_sim = lns_SN_200_mean + lns_SN_201_mean  # + lns_SN_202_mean
    labels_sim = [l.get_label() for l in handles_sim]
    plt.legend(handles_sim, labels_sim, loc='upper center', title='Simulations')

    ax_cal_degree.set_ylim(ymin=0)
    fname = fig_dir + 'stiffness_points_cal_degree'

    plt.savefig(fname)

    plt.show()
