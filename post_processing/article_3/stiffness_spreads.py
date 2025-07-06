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

matplotlib.style.use('axel_style_3')


def stiffness_spread_processing(sim_dir, sim_type, no_sims):
    simulation_directory = {}

    for i in range(1, no_sims + 1):
        simulation_directory[i] = MechanicalLoading(
            sim_dir + str(i) + '/swelling_electrode_mechanical_loading' + str(sim_type))

    stiffness_points_matrix = np.zeros((no_sims, simulation_directory[1].strain_points_total.size))
    strain_points = simulation_directory[1].strain_points_total
    for i in range(no_sims):
        stiffness_points_matrix[i, :] = simulation_directory[i + 1].stiffness_values_total

    stiffness_mean = np.mean(stiffness_points_matrix, axis=0)
    stiffness_std = np.std(stiffness_points_matrix, axis=0)
    return strain_points, stiffness_mean, stiffness_std


class MechanicalLoading:
    def __init__(self, sim_dir):
        self.sim_dir = sim_dir
        print(sim_dir)

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
        self.strain_points, self.stiffness_mean, self.stiffness_std = \
            stiffness_spread_processing(sim_dir, sim_type, no_sims)
        self.label = label


if __name__ == '__main__':

    # ==FIGURE SAVE DIR=================================================================================================
    fig_dir = 'C:/temp/figures/article_3/stiffness_spreads/'
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

    no_sims = 4
    SOC_0 = Simulation('/scratch/users/axlun/DEMsim/results/article_3/final_runs/', '', no_sims, '0')

    SOC_50 = Simulation('/scratch/users/axlun/DEMsim/results/article_3/final_runs/',
                        '_ss_1.03228_material_scaling', no_sims, '50')

    SOC_100 = Simulation('/scratch/users/axlun/DEMsim/results/article_3/final_runs/',
                         '_ss_1.06266_material_scaling', no_sims, '100')

    swelling = Simulation('/scratch/users/axlun/DEMsim/results/article_3/final_runs/',
                          '_ss_1.06266', no_sims, 'swelling')

    material_deg = Simulation('/scratch/users/axlun/DEMsim/results/article_3/final_runs/',
                              '_material_scaling', no_sims, 'material degradation')

    cycle_1 = Simulation('/scratch/users/axlun/DEMsim/results/article_3/final_runs/',
                         '_cycle_1', no_sims, '1')

    cycle_3 = Simulation('/scratch/users/axlun/DEMsim/results/article_3/final_runs/',
                         '_cycle_3', no_sims, '3')

    cycle_10 = Simulation('/scratch/users/axlun/DEMsim/results/article_3/final_runs/',
                          '_cycle_10', no_sims, '10')


    # =FIGURES==========================================================================================================
    # =FIGURE SOC=======================================================================================================
    fig_soc, ax_soc = plt.subplots()
    ax_soc.set_ylabel('Unloading stiffness [GPa]')
    ax_soc.set_xlabel('Strain [%]')
    ax_soc.set_ylim(ymin=0, ymax=5)
    ax_soc.xaxis.set_minor_locator(MultipleLocator(0.2))
    ax_soc.xaxis.set_major_locator(MultipleLocator(1))
    ax_soc.set_xlim(xmin=-2.2,xmax=2.2)
    # --SIMULATIONS----------------------------------------------------------------------------------------------------
    lns_soc_0_mean = ax_soc.plot(SOC_0.strain_points * 1E2, SOC_0.stiffness_mean * 1E-9,
                                 linestyle='dashed', marker='x', markersize=6, markeredgewidth=1, linewidth=1,
                                 color='C1', label=SOC_0.label)
    lns_soc_0_std = ax_soc.errorbar(SOC_0.strain_points * 1E2, SOC_0.stiffness_mean * 1E-9,
                                    yerr=SOC_0.stiffness_std * 1E-9, elinewidth=1, capthick=None, capsize=6,
                                    lw=0, marker='', markersize=0, markeredgewidth=1, color='C1',
                                    label='Standard deviation')

    lns_soc_50_mean = ax_soc.plot(SOC_50.strain_points * 1E2, SOC_50.stiffness_mean * 1E-9,
                                  linestyle='dashed', marker='x', markersize=6, markeredgewidth=1, linewidth=1,
                                  color='C2', label=SOC_50.label)
    lns_soc_50_std = ax_soc.errorbar(SOC_50.strain_points * 1E2, SOC_50.stiffness_mean * 1E-9,
                                     yerr=SOC_50.stiffness_std * 1E-9, elinewidth=1, capthick=None, capsize=6,
                                     lw=0, marker='', markersize=0, markeredgewidth=1, color='C2',
                                     label='Standard deviation')

    lns_soc_100_mean = ax_soc.plot(SOC_100.strain_points * 1E2, SOC_100.stiffness_mean * 1E-9,
                                   linestyle='dashed', marker='x', markersize=6, markeredgewidth=1, linewidth=1,
                                   color='C0', label=SOC_100.label)
    lns_soc_100_std = ax_soc.errorbar(SOC_100.strain_points * 1E2, SOC_100.stiffness_mean * 1E-9,
                                      yerr=SOC_100.stiffness_std * 1E-9, elinewidth=1, capthick=None, capsize=6,
                                      lw=0, marker='', markersize=0, markeredgewidth=1, color='C0',
                                      label='Standard deviation')

    handles_sim = lns_soc_100_mean + lns_soc_50_mean + lns_soc_0_mean
    labels_sim = [l.get_label() for l in handles_sim]
    plt.legend(handles_sim, labels_sim, loc='upper right', title='SOC [%]')
    fig_soc.tight_layout()
    fname = fig_dir + 'soc'
    plt.savefig(fname)
    plt.savefig(fname+'.svg')

    # =FIGURE SPLIT=======================================================================================================
    fig_split, ax_split = plt.subplots()
    ax_split.set_ylabel('Unloading stiffness [GPa]')
    ax_split.set_xlabel('Strain [%]')
    ax_split.set_ylim(ymin=0, ymax=5)
    ax_split.set_xlim(xmin=-2.2, xmax=2.2)
    ax_split.xaxis.set_minor_locator(MultipleLocator(0.2))
    ax_split.xaxis.set_major_locator(MultipleLocator(1))
    # --SIMULATIONS----------------------------------------------------------------------------------------------------
    lns_soc_0_mean = ax_split.plot(SOC_0.strain_points * 1E2, SOC_0.stiffness_mean * 1E-9,
                                   linestyle='dashed', marker='x', markersize=6, markeredgewidth=1, linewidth=1,
                                   color='C1', label='Reference')
    lns_soc_0_std = ax_split.errorbar(SOC_0.strain_points * 1E2, SOC_0.stiffness_mean * 1E-9,
                                      yerr=SOC_0.stiffness_std * 1E-9, elinewidth=1, capthick=None, capsize=6,
                                      lw=0, marker='', markersize=0, markeredgewidth=1, color='C1',
                                      label='Standard deviation')

    lns_swelling_mean = ax_split.plot(swelling.strain_points * 1E2, swelling.stiffness_mean * 1E-9,
                                      linestyle='dashed', marker='x', markersize=6, markeredgewidth=1, linewidth=1,
                                      color='C3', label='Swelling')
    lns_swelling_std = ax_split.errorbar(swelling.strain_points * 1E2, swelling.stiffness_mean * 1E-9,
                                         yerr=swelling.stiffness_std * 1E-9, elinewidth=1, capthick=None, capsize=6,
                                         lw=0, marker='', markersize=0, markeredgewidth=1, color='C3',
                                         label='Standard deviation')

    lns_mat_deg_mean = ax_split.plot(material_deg.strain_points * 1E2, material_deg.stiffness_mean * 1E-9,
                                     linestyle='dashed', marker='x', markersize=6, markeredgewidth=1, linewidth=1,
                                     color='C2', label='Material degradation')
    lns_mat_deg_std = ax_split.errorbar(material_deg.strain_points * 1E2, material_deg.stiffness_mean * 1E-9,
                                        yerr=material_deg.stiffness_std * 1E-9, elinewidth=1, capthick=None, capsize=6,
                                        lw=0, marker='', markersize=0, markeredgewidth=1, color='C2',
                                        label='Standard deviation')

    lns_soc_100_mean = ax_split.plot(SOC_100.strain_points * 1E2, SOC_100.stiffness_mean * 1E-9,
                                     linestyle='dashed', marker='x', markersize=6, markeredgewidth=1, linewidth=1,
                                     color='C0', label='Swelling + material degradation')
    lns_soc_100_std = ax_split.errorbar(SOC_100.strain_points * 1E2, SOC_100.stiffness_mean * 1E-9,
                                        yerr=SOC_100.stiffness_std * 1E-9, elinewidth=1, capthick=None, capsize=6,
                                        lw=0, marker='', markersize=0, markeredgewidth=1, color='C0',
                                        label='Standard deviation')

    handles_sim = lns_swelling_mean + lns_soc_100_mean + lns_mat_deg_mean + lns_soc_0_mean
    labels_sim = [l.get_label() for l in handles_sim]
    plt.legend(handles_sim, labels_sim, loc='lower left')
    fig_split.tight_layout()
    fname = fig_dir + 'split'
    plt.savefig(fname)
    plt.savefig(fname+'.svg')

    # =FIGURE CYCLING===================================================================================================
    fig_cycling, ax_cycling = plt.subplots()
    ax_cycling.set_ylabel('Unloading stiffness [GPa]')
    ax_cycling.set_xlabel('Strain [%]')
    ax_cycling.set_ylim(ymin=0, ymax=5)
    ax_cycling.xaxis.set_minor_locator(MultipleLocator(0.2))
    ax_cycling.xaxis.set_major_locator(MultipleLocator(1))
    ax_cycling.set_xlim(xmin=-2.2, xmax=2.2)
    # --SIMULATIONS----------------------------------------------------------------------------------------------------
    lns_cycle_0_mean = ax_cycling.plot(SOC_0.strain_points * 1E2, SOC_0.stiffness_mean * 1E-9,
                                       linestyle='dashed', marker='x', markersize=6, markeredgewidth=1, linewidth=1,
                                       color='C1', label=SOC_0.label)
    lns_cycle_0_std = ax_cycling.errorbar(SOC_0.strain_points * 1E2, SOC_0.stiffness_mean * 1E-9,
                                          yerr=SOC_0.stiffness_std * 1E-9, elinewidth=1, capthick=None, capsize=6,
                                          lw=0, marker='', markersize=0, markeredgewidth=1, color='C1',
                                          label='Standard deviation')

    lns_cycle_1_mean = ax_cycling.plot(cycle_1.strain_points * 1E2, cycle_1.stiffness_mean * 1E-9,
                                       linestyle='dashed', marker='x', markersize=6, markeredgewidth=1, linewidth=1,
                                       color='C4', label=cycle_1.label)
    lns_cycle_1_std = ax_cycling.errorbar(cycle_1.strain_points * 1E2, cycle_1.stiffness_mean * 1E-9,
                                          yerr=cycle_1.stiffness_std * 1E-9, elinewidth=1, capthick=None, capsize=6,
                                          lw=0, marker='', markersize=0, markeredgewidth=1, color='C4',
                                          label='Standard deviation')

    lns_cycle_3_mean = ax_cycling.plot(cycle_3.strain_points * 1E2, cycle_3.stiffness_mean * 1E-9,
                                       linestyle='dashed', marker='x', markersize=6, markeredgewidth=1, linewidth=1,
                                       color='C5', label=cycle_3.label)
    lns_cycle_3_std = ax_cycling.errorbar(cycle_3.strain_points * 1E2, cycle_3.stiffness_mean * 1E-9,
                                          yerr=cycle_3.stiffness_std * 1E-9, elinewidth=1, capthick=None, capsize=6,
                                          lw=0, marker='', markersize=0, markeredgewidth=1, color='C5',
                                          label='Standard deviation')

    lns_cycle_10_mean = ax_cycling.plot(cycle_10.strain_points * 1E2, cycle_10.stiffness_mean * 1E-9,
                                        linestyle='dashed', marker='x', markersize=6, markeredgewidth=1, linewidth=1,
                                        color='C6', label=cycle_10.label)
    lns_cycle_10_std = ax_cycling.errorbar(cycle_10.strain_points * 1E2, cycle_10.stiffness_mean * 1E-9,
                                           yerr=cycle_10.stiffness_std * 1E-9, elinewidth=1, capthick=None, capsize=6,
                                           lw=0, marker='', markersize=0, markeredgewidth=1, color='C6',
                                           label='Standard deviation')

    handles_sim = lns_cycle_10_mean + lns_cycle_3_mean + lns_cycle_1_mean + lns_cycle_0_mean
    labels_sim = [l.get_label() for l in handles_sim]
    plt.legend(handles_sim, labels_sim, loc='upper right', title='Charge cycles')

    fig_cycling.tight_layout()
    fname = fig_dir + 'cycling'
    plt.savefig(fname)
    plt.savefig(fname+'.svg')

    print('Showing plots')
    plt.show()
