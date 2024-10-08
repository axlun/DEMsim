import os
from Bertil_functions.Bertil_functions import *

import numpy as np
import re
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.style
import time
import pandas as pd
import shutil
import os

if __name__ == '__main__':

    # ==NATUAL PACKING =================================================================================================
    # simulation_directory = '/scratch/users/axlun/DEMsim/results/article_2/final_runs/particle_contact_model/SN_2_2/electrode_natural_packing_el_pl_binder_el_pl_particle/'

    # ==CALENDERING======================================================================================================
    # simulation_directory = '/scratch/users/axlun/DEMsim/results/article_2/final_runs/particle_contact_model/SN_2_5_3/electrode_calendering_el_pl_binder_el_pl_particle/'

    # ==MECHANICAL LOADING==============================================================================================
    # simulation_directory = '/scratch/users/axlun/DEMsim/results/article_2/final_runs/particle_contact_model/SN_1_2_4/electrode_mechanical_loading_el_pl_binder_el_pl_particle_compression/'

    # ==RELAXATION======================================================================================================
    # simulation_directory = '/scratch/users/axlun/DEMsim/results/article_2/final_runs_2/SN_2111/1/' \
    #                        'electrode_relaxation_el_pl_binder_el_pl_particle_compression'

    # ==RELAXATION======================================================================================================
    # simulation_directory = '/scratch/users/axlun/DEMsim/results/article_2/viscoelastic_testing/SN_1/electrode_relaxation_el_pl_binder_el_pl_particle_compression'


    # simulation_directory = '/scratch/users/axlun/DEMsim/results/article_3/swelling/SN_2/swelling_electrode_calendering/'
    # simulation_directory = '/scratch/users/axlun/DEMsim/results/article_3/swelling/SN_2/electrode_swelling/'
    # simulation_directory = '/scratch/users/axlun/DEMsim/results/article_3/swelling/SN_2/swelling_electrode_mechanical_loading_tension/'
    # simulation_directory = '/scratch/users/axlun/DEMsim/results/article_3/swelling/SN_2/swelling_electrode_mechanical_loading_compression/'

    simulation_directory = '/scratch/users/axlun/DEMsim/results/article_3/swelling_material_scaling/SN_1/electrode_cycling_1/'

    # ==FIGURE DIR======================================================================================================
    fig_dir = 'C:/temp/figures/Bertil_sum_normal_tangential_force/'
    try:
        shutil.rmtree(fig_dir)
        os.mkdir(fig_dir)
    except FileNotFoundError:
        print('Directory does not exist')
        os.mkdir(fig_dir)
    except:
        print('Directory could not be removed')
        quit()

    time_steps, normal_force, compressive_normal_force_vec, tensile_normal_force_vec, binder_force_vec, \
    compressive_binder_force_vec, tensile_binder_force_vec, particle_force_vec, tangential_force = \
        sum_normal_tangential_force_bertil(simulation_directory)

    # ==PLOTTING PARAMETERS=============================================================================================
    plt.style.use('axel_style')

    # ==FIG 1 NORMAL AND TANGENTIAL FORCES==============================================================================
    figure_normal_tangential_force_time, ax_normal_tangential_force_time = plt.subplots()
    lns_normal_force_time = ax_normal_tangential_force_time.plot(time_steps[:], normal_force, 'r', linewidth=3,
                                                                 label=r'$F_n$')
    lns_tangential_force_time = ax_normal_tangential_force_time.plot(time_steps[:], tangential_force, 'g', linewidth=3,
                                                                     label=r'$F_t$')
    ax_normal_tangential_force_time.set_xlabel('Time [s]')
    ax_normal_tangential_force_time.set_ylabel('Total force [N]')
    ax_normal_tangential_force_time.set_ylim(ymin=0)#,ymax=10)
    ax_normal_tangential_force_time.legend()
    figure_normal_tangential_force_time.tight_layout()

    fname = fig_dir + 'normal_tangential_force'
    plt.savefig(fname)


    # ==FIG 2 NORMAL FORCE IN COMPRESSION AND TENSION===================================================================
    figure_normal_compression_tension_force_time, ax_normal_compression_tension_force_time = plt.subplots()
    lns_normal_force = ax_normal_compression_tension_force_time.plot(time_steps[:], normal_force, 'r', linewidth=3,
                                                                     label=r'$F_n$')
    lns_normal_compression_force = ax_normal_compression_tension_force_time.plot(time_steps[:],
                                                                                 compressive_normal_force_vec, 'g',
                                                                                 linewidth=3,
                                                                                 label=r'$F_{n,comp}$')
    lns_normal_tensile_force = ax_normal_compression_tension_force_time.plot(time_steps[:],
                                                                             tensile_normal_force_vec, 'b',
                                                                             linewidth=3,
                                                                             label=r'$F_{n,ten}$')

    ax_normal_compression_tension_force_time.set_xlabel('Time [s]')
    ax_normal_compression_tension_force_time.set_ylabel('Total force [N]')
    ax_normal_compression_tension_force_time.set_ylim(ymin=0)
    ax_normal_compression_tension_force_time.legend()
    figure_normal_compression_tension_force_time.tight_layout()

    fname = fig_dir + 'tensile_compressive_normal_force'
    plt.savefig(fname)

    # ==FIG 3 BINDER FORCES=============================================================================================
    figure_binder_forces, ax_binder_forces = plt.subplots()
    lns_binder_force_time = ax_binder_forces.plot(time_steps, binder_force_vec, 'r', linewidth=3, label=r'$F_b$')
    lns_compressive_binder_force_time = ax_binder_forces.plot(time_steps, compressive_binder_force_vec, 'g',
                                                              linewidth=3, label=r'$F_{b,comp}$')
    lns_tensile_binder_force_time = ax_binder_forces.plot(time_steps, tensile_binder_force_vec, 'b', linewidth=3,
                                                          label=r'$F_{b,ten}$')
    ax_binder_forces.set_xlabel('Time [s]')
    ax_binder_forces.set_ylabel('Total force [N]')
    ax_binder_forces.set_ylim(ymin=0)
    ax_binder_forces.legend()
    figure_binder_forces.tight_layout()

    fname = fig_dir + 'binder_normal_forces'
    plt.savefig(fname)

    # ==FIG 4 TOTAL BINDER AND PARTICLE FORCES==========================================================================
    figure_total_binder_particle_forces, ax_total_binder_particle_forces = plt.subplots()
    lns_normal_force_time2 = ax_total_binder_particle_forces.plot(time_steps, normal_force, 'r', linewidth=3,
                                                                  label=r'$F_n$')
    lns_binder_force_time2 = ax_total_binder_particle_forces.plot(time_steps, binder_force_vec, 'g',
                                                              linewidth=3, label=r'$F_b$')
    lns_particle_force_time = ax_total_binder_particle_forces.plot(time_steps, particle_force_vec, 'b', linewidth=3,
                                                          label=r'$F_p$')
    ax_total_binder_particle_forces.set_xlabel('Time [s]')
    ax_total_binder_particle_forces.set_ylabel('Total force [N]')
    ax_total_binder_particle_forces.set_ylim(ymin=0)
    ax_total_binder_particle_forces.legend()
    figure_total_binder_particle_forces.tight_layout()

    fname = fig_dir + 'particle_binder_normal_forces'
    plt.savefig(fname)


    print("Plotting results")
    plt.show()
