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
matplotlib.style.use('axel_style')


class Calendering:
    def __init__(self, sim_dir, label):
        self.sim_dir = sim_dir

        self.time, self.calendering_surface_pressure, self.bottom_surface_pressure, self.calendering_surface_position, \
        self.sxx, self.syy, self.szz, self.tau_xy, self.tau_xz, self.tau_yx, self.tau_yz, self.tau_zx, self.tau_zy, \
        self.ke = calendering_plot_processing(sim_dir)
        self.label = label


class MechanicalLoading:
    def __init__(self, sim_dir, label):
        self.sim_dir = sim_dir

        self.time_tension, self.linear_strain_tension, self.sxx_tension, self.syy_tension, self.szz_tension,\
            self.tau_xy_tension, self.tau_xz_tension, self.tau_yx_tension, self.tau_yz_tension, self.tau_zx_tension,\
            self.tau_zy_tension, self.time_compression, self.linear_strain_compression, self.sxx_compression,\
            self.syy_compression, self.szz_compression, self.tau_xy_compression, self.tau_xz, self.tau_yx,\
            self.tau_yz_compression, self.tau_zx_compression, self.tau_zy_compression= mech_plot_prop(sim_dir)

        self.strain_points_total, self.stiffness_values_total = stiffness_func(self.sxx_tension,
                                                                              self.linear_strain_tension,
                                                                              self.sxx_compression,
                                                                              self.linear_strain_compression)
        self.label = label


class Simulation:
    def __init__(self, sim_dir, sim_type,label):
        self.sim_dir = sim_dir
        self.sim_type = sim_type
        self.calendering = Calendering(sim_dir + '/electrode_calendering_' + sim_type,label)
        self.mechanical_loading = MechanicalLoading(sim_dir + '/electrode_mechanical_loading_'+ sim_type,label)


if __name__=='__main__':

    # simulation_directory = '/scratch/users/axlun/DEMsim/results/article_2/final_runs/ref_sim/SN_1_2/'
    # simulation_1 = Simulation(simulation_directory, 'hertz')

    simulation_directory_H_SN_1 = '/scratch/users/axlun/DEMsim/results/article_2/final_runs/ref_sim/SN_1/'
    H_SN_1= Simulation(simulation_directory_H_SN_1, 'hertz', 'SN_1')
    # ==FIGURE SAVE DIR=================================================================================================
    fig_dir = 'C:/temp/figures/article_2/calendering_and_stiffness/'
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


    # ==EXPERIMENTAL DATA===================================================================================================
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


    # =FIGURE 1=========================================================================================================
    fig_calendering, ax_calendering = plt.subplots()
    ax_calendering.plot(H_SN_1.calendering.calendering_surface_position[:] * 1E2,
                        H_SN_1.calendering.calendering_surface_pressure[:] * 1E-6,
                        label=H_SN_1.calendering.label)
    ax_calendering.set_xlim(xmin=104.8,xmax=120)
    ax_calendering.xaxis.set_major_locator(MultipleLocator(5))
    ax_calendering.set_ylim(ymin=0)
    ax_calendering.set_ylabel('Calendering surface pressure [MPa]')
    ax_calendering.set_xlabel('Calendering surface hetight [µm]')
    fig_calendering.tight_layout()
    ax_calendering.legend(loc='best')
    fname=fig_dir + 'calendering_surface_pressure_to_position'
    plt.savefig(fname)

    # =FIGURE 2=========================================================================================================
    fig_stiffness, ax_stiffness = plt.subplots()
    ax_stiffness.set_ylabel('Unloading stiffness [GPa]')
    ax_stiffness.set_xlabel('Strain [%]')
    # ax_stiff_El_Pl.set_title('Unloading stiffness of electrode layer')
    # --SIMULATIONS----------------------------------------------------------------------------------------------------
    lns_stiff= ax_stiffness.plot(
        H_SN_1.mechanical_loading.strain_points_total * 1E2,
        H_SN_1.mechanical_loading.stiffness_values_total * 1E-9,
        linestyle='dashed', marker='x', markersize=12, markeredgewidth=3, linewidth=3,
        label=H_SN_1.mechanical_loading.label)

    #lns_stiff_El_Pl = ax_stiffness.plot(strain_points_total_SN_run_1_El_Pl * 100,
    #                                    stiffness_values_total_SN_run_1_El_Pl * 1E-9, linestyle='dashed',
    #                                    marker='o', fillstyle='none', markersize=12, markeredgewidth=3, linewidth=3,
    #                                    label=r'$El-Pl$')

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
    handles_exp = lns_stiff_exp_eps_dot_01_El_Pl + lns_stiff_exp_eps_dot_05_El_Pl + lns_stiff_exp_eps_dot_10_El_Pl + lns_stiff_exp_eps_dot_100_El_Pl + lns_stiff_exp_eps_dot_300_El_Pl
    handles_exp.reverse()
    labels_exp = [l.get_label() for l in handles_exp]
    first_legend_El_Pl = plt.legend(handles_exp, labels_exp, loc='upper right', title='Experimental results')

    ax_stiffness.add_artist(first_legend_El_Pl)

    handles_sim= lns_stiff #+ lns_stiff_El_Pl
    labels_sim = [l.get_label() for l in handles_sim]
    plt.legend(handles_sim, labels_sim, loc='upper center', title='Simulations')

    ax_stiffness.set_ylim(ymin=0)
    fname = fig_dir + 'stiffness_points'
    plt.savefig(fname)

    plt.show()