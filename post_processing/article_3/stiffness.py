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

""""
def stiffness_processing(sim_dir):
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
"""


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
    def __init__(self, sim_dir, label='template'):
        self.sim_dir = sim_dir
        self.stress_stiffness = MechanicalLoading(self.sim_dir)
        self.strain_points = self.stress_stiffness.strain_points_total
        self.stiffness_points = self.stress_stiffness.stiffness_values_total
        # self.strain_points, self.stiffness_points = stiffness_processing(sim_dir)
        self.label = label


if __name__ == '__main__':

    # ==FIGURE SAVE DIR=================================================================================================
    fig_dir = 'C:/temp/figures/stiffness/'
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

    SOC_0 = Simulation('/scratch/users/axlun/DEMsim/results/article_3/charge_cycling/SN_1/'
                       'swelling_electrode_mechanical_loading', '0 %')

    SOC_50 = Simulation('/scratch/users/axlun/DEMsim/results/article_3/charge_cycling/SN_1/'
                        'swelling_electrode_mechanical_loading_ss_1.03228_material_scaling', '50 %')

    SOC_100 = Simulation('/scratch/users/axlun/DEMsim/results/article_3/charge_cycling/SN_1/'
                         'swelling_electrode_mechanical_loading_ss_1.06266_material_scaling', '100 %')

    cycle_1 = Simulation('/scratch/users/axlun/DEMsim/results/article_3/charge_cycling/SN_1/'
                         'swelling_electrode_mechanical_loading_cycle_1', '1')

    cycle_3 = Simulation('/scratch/users/axlun/DEMsim/results/article_3/charge_cycling/SN_1/'
                         'swelling_electrode_mechanical_loading_cycle_3', '3')

    cycle_10 = Simulation('/scratch/users/axlun/DEMsim/results/article_3/charge_cycling/SN_1/'
                          'swelling_electrode_mechanical_loading_cycle_10', '10')

    swelling = Simulation('/scratch/users/axlun/DEMsim/results/article_3/charge_cycling/SN_1/'
                          'swelling_electrode_mechanical_loading_ss_1.06266', 'swelling')

    material_scaling = Simulation('/scratch/users/axlun/DEMsim/results/article_3/charge_cycling/SN_1/'
                                  'swelling_electrode_mechanical_loading_material_scaling', 'material degradation')

    material_scaling_05 = Simulation('/scratch/users/axlun/DEMsim/results/article_3/charge_cycling/SN_1/'
                                     'swelling_electrode_mechanical_loading_material_scaling_05',
                                     'material degradation 50 %')

    # =FIGURES==========================================================================================================

    # =FIGURE SOC ======================================================================================================
    fig_Sep, ax_Sep = plt.subplots()
    ax_Sep.set_ylabel('Unloading stiffness [GPa]')
    ax_Sep.set_xlabel('Strain [%]')
    # --SIMULATIONS----------------------------------------------------------------------------------------------------
    lns_SOC_0 = ax_Sep.plot(SOC_0.strain_points * 1E2, SOC_0.stiffness_points * 1E-9,
                            linestyle='dashed', marker='x', markersize=12, markeredgewidth=3, linewidth=3,
                            color='C0', label='unchanged')

    lns_swelling = ax_Sep.plot(swelling.strain_points * 1E2, swelling.stiffness_points * 1E-9,
                                 linestyle='dashed', marker='x', markersize=12, markeredgewidth=3, linewidth=3,
                                 color='C1', label=swelling.label)

    lns_material_scaling = ax_Sep.plot(material_scaling.strain_points * 1E2, material_scaling.stiffness_points * 1E-9,
                                       linestyle='dashed', marker='x', markersize=12, markeredgewidth=3, linewidth=3,
                                       color='C2', label=material_scaling.label)

    lns_material_scaling_05 = ax_Sep.plot(material_scaling_05.strain_points * 1E2, material_scaling_05.stiffness_points * 1E-9,
                                       linestyle='dashed', marker='x', markersize=12, markeredgewidth=3, linewidth=3,
                                          color='C4', label=material_scaling_05.label)

    lns_SOC_100 = ax_Sep.plot(SOC_100.strain_points * 1E2, SOC_100.stiffness_points * 1E-9,
                              linestyle='dashed', marker='x', markersize=12, markeredgewidth=3, linewidth=3,
                              color = 'C3', label = 'swelling + material degradation')

    handles_sim = lns_SOC_100 + lns_swelling + lns_material_scaling + lns_material_scaling_05 + lns_SOC_0
    labels_sim = [l.get_label() for l in handles_sim]
    plt.legend(handles_sim, labels_sim, loc='lower left')

    ax_Sep.set_ylim(ymin=0, ymax=4.5)
    fname = fig_dir + 'Sep'
    plt.savefig(fname)

    # =FIGURE SOC ======================================================================================================
    fig_SOC, ax_SOC = plt.subplots()
    ax_SOC.set_ylabel('Unloading stiffness [GPa]')
    ax_SOC.set_xlabel('Strain [%]')
    # --SIMULATIONS----------------------------------------------------------------------------------------------------
    lns_SOC_0 = ax_SOC.plot(SOC_0.strain_points * 1E2, SOC_0.stiffness_points * 1E-9,
                            linestyle='dashed', marker='x', markersize=12, markeredgewidth=3, linewidth=3,
                            color='C0', label=SOC_0.label)
    lns_SOC_50 = ax_SOC.plot(SOC_50.strain_points * 1E2, SOC_50.stiffness_points * 1E-9,
                             linestyle='dashed', marker='x', markersize=12, markeredgewidth=3, linewidth=3,
                             color='C1', label=SOC_50.label)

    lns_SOC_100 = ax_SOC.plot(SOC_100.strain_points * 1E2, SOC_100.stiffness_points * 1E-9,
                              linestyle='dashed', marker='x', markersize=12, markeredgewidth=3, linewidth=3,
                              color='C3', label=SOC_100.label)

    handles_sim = lns_SOC_100 + lns_SOC_50 + lns_SOC_0
    labels_sim = [l.get_label() for l in handles_sim]
    plt.legend(handles_sim, labels_sim, loc='lower left', title='SOC')  # , title='Simulations')

    ax_SOC.set_ylim(ymin=0, ymax=4)
    fname = fig_dir + 'SOC'
    plt.savefig(fname)

    # =FIGURE CYCLES====================================================================================================
    fig_cycle, ax_cycle = plt.subplots()
    ax_cycle.set_ylabel('Unloading stiffness [GPa]')
    ax_cycle.set_xlabel('Strain [%]')
    # --SIMULATIONS----------------------------------------------------------------------------------------------------
    lns_cycle_0 = ax_cycle.plot(SOC_0.strain_points * 1E2, SOC_0.stiffness_points * 1E-9,
                                linestyle='dashed', marker='x', markersize=12, markeredgewidth=3, linewidth=3,
                                color='C0', label='0')

    lns_cycle_1 = ax_cycle.plot(cycle_1.strain_points * 1E2, cycle_1.stiffness_points * 1E-9,
                                linestyle='dashed', marker='x', markersize=12, markeredgewidth=3, linewidth=3,
                                color='C1', label=cycle_1.label)

    lns_cycle_3 = ax_cycle.plot(cycle_3.strain_points * 1E2, cycle_3.stiffness_points * 1E-9,
                                linestyle='dashed', marker='x', markersize=12, markeredgewidth=3, linewidth=3,
                                color='C2', label=cycle_3.label)

    lns_cycle_10 = ax_cycle.plot(cycle_10.strain_points * 1E2, cycle_10.stiffness_points * 1E-9,
                                 linestyle='dashed', marker='x', markersize=12, markeredgewidth=3, linewidth=3,
                                 color='C3', label=cycle_10.label)

    handles_sim = lns_cycle_10 + lns_cycle_3 + lns_cycle_1 + lns_cycle_0
    labels_sim = [l.get_label() for l in handles_sim]
    plt.legend(handles_sim, labels_sim, loc='lower left', title='Cycles')  # , title='Simulations')

    ax_cycle.set_ylim(ymin=0, ymax=4)
    fname = fig_dir + 'cycles'
    plt.savefig(fname)

    plt.show()
