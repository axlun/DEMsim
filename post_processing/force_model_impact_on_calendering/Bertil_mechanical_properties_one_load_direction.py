from Bertil_calendering_pressure import local_data_gatherer, bertil_data_gatherer
from Bertil_mechanical_properties import stress_and_linear_strain_finder,stiffness_finder
import numpy as np

import matplotlib.pyplot as plt
import matplotlib
import os
import shutil

if __name__ == '__main__':

    # simulation_directory = '/scratch/users/axlun/DEMsim/results/electrode_mechanical_loading_hertz/SN_hertz_5000p_btr_5_brr_05_dt_5e1_MS_1e2_SR_2e-3_tension'

    # simulation_directory = '/scratch/users/axlun/DEMsim/results/electrode_mechanical_loading_hertz/SN_ref_run_1_5000p_btr_5_brr_15_dt_5e1_MS_1e2_SR_2e-3_no_new_binder_run_3_tension'
    # simulation_directory = '/scratch/users/axlun/DEMsim/results/electrode_mechanical_loading_hertz/SN_ref_run_1_5000p_btr_5_brr_15_dt_5e1_MS_1e2_SR_2e-3_rot_compression'
    # simulation_directory = '/scratch/users/axlun/DEMsim/results/final_runs/SN_run_1_br_175/electrode_mechanical_loading_hertz_compression'
    simulation_directory = '/scratch/users/axlun/DEMsim/results/final_runs/SN_run_1_q_0_rmin_3_rmax_10/electrode_mechanical_loading_hertz_tension'



    stiffness_at_points_flag = 1

    fig_dir = 'C:/temp/figures/Bertil_mechanical_properties_one_load_direction/'
    try:
        shutil.rmtree(fig_dir)
        os.mkdir(fig_dir)
    except FileNotFoundError:
        print('Directory does not exist')
        os.mkdir(fig_dir)
    except:
        print('Directory could not be removed')
        quit()

    # ==PLOT PARAMETERS=================================================================================================
    plt.style.use('axel_style')

    # ==EXPERIMENTAL DATA===============================================================================================

    exp_strain_points_compression = [-1.00, -1.10, -1.23, -1.41, -1.65]
    Modulus_eps_dot_01_compression = [1.20, 1.04, 0.96, 1.19, 1.46]
    Modulus_eps_dot_05_compression = [1.43, 0.86, 0.97, 1.39, 1.79]
    Modulus_eps_dot_10_compression = [1.51, 1.79, 1.99, 2.09, 2.74]
    Modulus_eps_dot_100_compression = [1.55, 1.72, 1.70, 1.81, 2.28]
    Modulus_eps_dot_300_compression = [1.99, 1.52, 2.21, 2.80, 2.33]

    exp_strain_points_tension = [1.00, 1.10, 1.23, 1.41, 1.65]
    Modulus_eps_dot_01_tension = [0.95, 0.50, 0.89, 0.61, 0.57]
    Modulus_eps_dot_05_tension = [0.78, 0.39, 1.17, 0.55, 0.46]
    Modulus_eps_dot_10_tension = [0.84, 0.59, 1.09, 0.59, 0.52]
    Modulus_eps_dot_100_tension = [0.76, 0.55, 0.91, 0.63, 0.67]
    Modulus_eps_dot_300_tension = [1.06, 0.55, 0.90, 0.95, 0.85]

    tension = False
    compression = False
    if 'tension' in simulation_directory:
        tension = True
    else:
        compression = True
    periodic_BC_data_tension, force_fabric_tensor_data_tension, surface_position_data_tension = [],[],[]
    if simulation_directory.startswith("/scratch"):
        force_data_tension, surface_force_index_tension, surface_position_index_tension, surface_position_data_tension, periodic_BC_data_tension, force_fabric_tensor_data_tension, kinetic_energy_data_tension = bertil_data_gatherer(
            simulation_directory)
    else:
        print("Error with simulation directory")
    time, linear_strain, sxx, syy, szz, tau_xy, tau_xz, tau_yx, tau_yz, tau_zx, tau_zy = stress_and_linear_strain_finder(periodic_BC_data_tension,
                                                                         force_fabric_tensor_data_tension,
                                                                         surface_position_data_tension)

    # ==FIG 1 STRESS XX YY ZZ STRAIN====================================================================================
    fig_tension,ax_tension = plt.subplots()
    ax_tension.set_ylabel('Stress [MPa]')
    ax_tension.set_xlabel('Strain [%]')
    lns_tension_sig_xx = ax_tension.plot(linear_strain[:] * 100, sxx[:] / 1e6, label=r'$\sigma_{xx}$')
    lns__tension_sig_yy = ax_tension.plot(linear_strain[:] * 100, syy[:] / 1e6, label=r'$\sigma_{yy}$')
    lns_tension_sig_zz = ax_tension.plot(linear_strain[:] * 100, szz[:] / 1e6, label=r'$\sigma_{zz}$')

    lns_tension_tau_xy = ax_tension.plot(linear_strain[:] * 100, tau_xy[:] / 1e6, label=r'$\tau_{xy}$')
    lns_tension_tau_xz = ax_tension.plot(linear_strain[:] * 100, tau_xz[:] / 1e6, label=r'$\tau_{xz}$')
    lns_tension_tau_yx = ax_tension.plot(linear_strain[:] * 100, tau_yx[:] / 1e6, label=r'$\tau_{yx}$')
    lns_tension_tau_yz = ax_tension.plot(linear_strain[:] * 100, tau_yz[:] / 1e6, label=r'$\tau_{yz}$')
    lns_tension_tau_zx = ax_tension.plot(linear_strain[:] * 100, tau_zx[:] / 1e6, label=r'$\tau_{zx}$')
    lns_tension_tau_zy = ax_tension.plot(linear_strain[:] * 100, tau_zy[:] / 1e6, label=r'$\tau_{zy}$')

    ax_tension.set_title('Macroscopic stress')
    plt.legend(loc='best')
    fig_tension.tight_layout()

    fname = fig_dir + 'stress_strain'
    plt.savefig(fname)

    # ==FIG 2 STRESS STRAIN TIME =======================================================================================
    fig_tension_time, ax_tension_time = plt.subplots()
    ax_tension_time.set_ylabel("Stress [MPa]")
    ax_tension_time.set_xlabel("time [s]")
    lns_tension_xx_time = ax_tension_time.plot(time, sxx[:] / 1e6, label=r'$\sigma_{xx}$')
    lns_tension_yy_time = ax_tension_time.plot(time, syy[:] / 1e6,
                                               label=r'$\sigma_{yy}$')
    lns_tension_zz_time = ax_tension_time.plot(time, szz[:] / 1e6,
                                               label=r'$\sigma_{zz}$')
    # ax_tension_time.set_ylim(ymin=-10,ymax=10)

    ax_tension_time_2 = ax_tension_time.twinx()
    ##    ax2.set_xlabel("Calendering surface position [m]")
    #    ax2.set_ylabel("Calendering surface position [m]")
    lns_tension_time_2 = ax_tension_time_2.plot(time, linear_strain * 100,'c', label=r'$\varepsilon_{xx}$')
    ax_tension_time_2.set_ylabel("Strain [%]")
    # lns3 = ax2.plot(time[max_index], linear_strain[max_index],'b',label='Position',linestyle="None",marker='x',markeredgewidth = 3)
    # lns4 = ax2.plot(time[min_index], linear_strain[min_index], 'b', label='Position', linestyle="None", marker='x',markeredgewidth = 3)

    lns_tension = lns_tension_xx_time +lns_tension_yy_time+lns_tension_zz_time+ lns_tension_time_2
    labs_tension = [l.get_label() for l in lns_tension]
    ax_tension_time.legend(lns_tension, labs_tension, loc=0)
    ax_tension_time.set_title('Macroscopic stress and strain')
    fig_tension_time.tight_layout()

    fname = fig_dir + 'stress_strain_time'
    plt.savefig(fname)


    # ==FIG 3 STIFFNESS AT STRAIN POINTS================================================================================

    if stiffness_at_points_flag == 1:
        strain_points_tension, stiffness_values_tension = stiffness_finder(sxx, linear_strain)

        strain_points_total = strain_points_tension
        stiffness_values_total = stiffness_values_tension

        fig_stiff, ax_stiff = plt.subplots()
        ax_stiff.set_ylabel('Stiffness [GPa]')
        ax_stiff.set_xlabel('Strain [%]')
        ax_stiff.set_title('Stiffness of electrode layer')
        lns_stiff = ax_stiff.plot(strain_points_total * 100, stiffness_values_total * 1E-9, color='red', linestyle='dashed',
                                  marker='x',markerfacecolor='blue', markersize=12, markeredgewidth=3,
                                  linewidth=3, label='Simulation')

        if tension:
            lns_stiff_exp_tension_eps_dot_01 = ax_stiff.plot(exp_strain_points_tension,Modulus_eps_dot_01_tension,
                                                             marker='x',markersize=12, markeredgewidth=3,
                                                             linewidth=0,label='0.1 mm/min')

            lns_stiff_exp_tension_eps_dot_05 = ax_stiff.plot(exp_strain_points_tension,Modulus_eps_dot_05_tension,
                                                             marker='x',markersize=12, markeredgewidth=3,
                                                             linewidth=0,label='0.5 mm/min')
            lns_stiff_exp_tension_eps_dot_10 = ax_stiff.plot(exp_strain_points_tension, Modulus_eps_dot_10_tension,
                                                             marker='x', markersize=12, markeredgewidth=3,
                                                             linewidth=0, label='1.0 mm/min')
            lns_stiff_exp_tension_eps_dot_100 = ax_stiff.plot(exp_strain_points_tension, Modulus_eps_dot_100_tension,
                                                             marker='x', markersize=12, markeredgewidth=3,
                                                             linewidth=0, label='10 mm/min')
            lns_stiff_exp_tension_eps_dot_300 = ax_stiff.plot(exp_strain_points_tension, Modulus_eps_dot_300_tension,
                                                             marker='x', markersize=12, markeredgewidth=3,
                                                             linewidth=0, label='30 mm/min')

        if compression:
            lns_stiff_exp_compression_eps_dot_01 = ax_stiff.plot(exp_strain_points_compression, Modulus_eps_dot_01_compression,
                                                             marker='x', markersize=12, markeredgewidth=3,
                                                             linewidth=0, label='0.1 mm/min')

            lns_stiff_exp_compression_eps_dot_05 = ax_stiff.plot(exp_strain_points_compression, Modulus_eps_dot_05_compression,
                                                             marker='x', markersize=12, markeredgewidth=3,
                                                             linewidth=0, label='0.5 mm/min')
            lns_stiff_exp_compression_eps_dot_10 = ax_stiff.plot(exp_strain_points_compression, Modulus_eps_dot_10_compression,
                                                             marker='x', markersize=12, markeredgewidth=3,
                                                             linewidth=0, label='1.0 mm/min')
            lns_stiff_exp_compression_eps_dot_100 = ax_stiff.plot(exp_strain_points_compression, Modulus_eps_dot_100_compression,
                                                              marker='x', markersize=12, markeredgewidth=3,
                                                              linewidth=0, label='10 mm/min')
            lns_stiff_exp_compression_eps_dot_300 = ax_stiff.plot(exp_strain_points_compression, Modulus_eps_dot_300_compression,
                                                              marker='x', markersize=12, markeredgewidth=3,
                                                              linewidth=0, label='30 mm/min')


        ax_stiff.set_ylim(ymin=0)
        fig_stiff.tight_layout()

        plt.legend(loc='best')
        fname = fig_dir + 'stiffness_strain'
        plt.savefig(fname)
    # ==FIG 4 KINETIC ENERGY============================================================================================
    fig_KE,ax_KE = plt.subplots()
    ax_KE.set_ylabel('Kinetic energy [J]')
    ax_KE.set_xlabel('Time [s]')
    lns_KE = ax_KE.plot(kinetic_energy_data_tension[:, -1].astype(float),
                        kinetic_energy_data_tension[:, 0].astype(float), label=r'KE')
    ax_KE.set_title('Kinetic energy')
    # plt.legend(loc='best')
    fig_KE.tight_layout()
    fname = fig_dir + 'KE'
    plt.savefig(fname)

    # ==SHOWIN PLOT=====================================================================================================
    plt.show()
