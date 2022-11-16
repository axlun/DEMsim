from Bertil_calendering_pressure import local_data_gatherer, bertil_data_gatherer

import numpy as np
import matplotlib.pyplot as plt
import shutil
import os
import matplotlib
# matplotlib.style.use('classic')


def stress_and_linear_strain_finder(periodic_BC_data, force_fabric_tensor_data, surface_position_data):
    time = periodic_BC_data[:, 0]

    x_side_length = periodic_BC_data[:, 2] - periodic_BC_data[:, 1]
    x_side_length_0 = x_side_length[0]
    y_side_length = periodic_BC_data[:, 4] - periodic_BC_data[:, 3]

    #t0 = 1.11  # Chould be chaged later, defined in another way?======================================================
    t0 = float(surface_position_data[0,32])-1
    vol = x_side_length * y_side_length * t0
    sxx = -force_fabric_tensor_data[:, 1] / vol
    syy = -force_fabric_tensor_data[:, 5] / vol
    szz = -force_fabric_tensor_data[:, 9] / vol
    linear_strain = (x_side_length[:] - x_side_length_0) / x_side_length_0

    return time, linear_strain, sxx, syy, szz


def stiffness_finder(stress_vec, strain_vec):
    count = 0
    n = len(strain_vec)
    max_index = []
    min_index = []
    # Tension
    if strain_vec[1] > 0:
        for i in range(1, n - 1):
            # maximum values
            if (strain_vec[i] > strain_vec[i + 1] and strain_vec[i] >= strain_vec[i - 1]):
                max_index.append(i)
                count += 1
            # Minimum values
            if (strain_vec[i] <= strain_vec[i + 1] and strain_vec[i] < strain_vec[i - 1]):
                min_index.append(i)
                count += 1
            if count >= 10:
                break
    # Compression
    if strain_vec[1] < 0:
        for i in range(1, n - 1):
            # Maximum values
            if (strain_vec[i] >= strain_vec[i + 1] and strain_vec[i] > strain_vec[i - 1]):
                max_index.append(i)
                count += 1
            # Minimum values
            if (strain_vec[i] < strain_vec[i + 1] and strain_vec[i] <= strain_vec[i - 1]):
                min_index.append(i)
                count += 1
            if count >= 10:
                break
    stiffness_values = np.array([])
    for i in range(5):
        # ==================STRAIN VÄRDEN LÄGGS INTE TILL I ARRAYN...===========
        stiffness_values = np.append(stiffness_values, (stress_vec[max_index[i]] - stress_vec[min_index[i]]) / (strain_vec[max_index[i]] - strain_vec[min_index[i]]))
    if strain_vec[1] < 0:
        strain_points = strain_vec[min_index]
    if strain_vec[1] > 0:
        strain_points = strain_vec[max_index]

    return strain_points, stiffness_values  # returns strain in [-] and stiffness in Pa

if __name__ == '__main__':

    simulation_directory = '/scratch/users/axlun/DEMsim/results/electrode_mechanical_loading_hertz/SN_ref_run_1_5000p_btr_5_brr_15_dt_5e1_MS_1e2_SR_2e-3'

    stiffness_at_points_flag = 1

    fig_dir = 'C:/temp/figures/Bertil_mechanical_properties/'
    try:
        shutil.rmtree(fig_dir)
        os.mkdir(fig_dir)
    except FileNotFoundError:
        print('Directory does not exist')
        os.mkdir(fig_dir)
    except:
        print('Directory could not be removed')
        quit()
# ==PLOT PARAMETERS=====================================================================================================
    plt.rcParams['figure.figsize'] = (12,9)
    plt.rcParams['lines.linewidth'] = 2
    plt.rcParams['axes.labelweight'] = 'bold'
    plt.rcParams['axes.titleweight'] = 'bold'
    plt.rcParams['font.weight'] = 'bold'
    plt.rcParams['font.size'] = 20
    plt.rcParams["font.family"] = "Times New Roman"
    plt.rcParams['mathtext.default'] = 'regular'

# ==EXPERIMENTAL DATA===================================================================================================
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

    exp_strain_points = exp_strain_points_compression[::-1] + exp_strain_points_tension
    Modulus_eps_dot_01 = Modulus_eps_dot_01_compression[::-1] + Modulus_eps_dot_01_tension
    Modulus_eps_dot_05 = Modulus_eps_dot_05_compression[::-1] + Modulus_eps_dot_05_tension
    Modulus_eps_dot_10 = Modulus_eps_dot_10_compression[::-1] + Modulus_eps_dot_10_tension
    Modulus_eps_dot_100 = Modulus_eps_dot_100_compression[::-1] + Modulus_eps_dot_100_tension
    Modulus_eps_dot_300 = Modulus_eps_dot_300_compression[::-1] + Modulus_eps_dot_300_tension

    if simulation_directory.startswith("/scratch"):
        force_data_compression, surface_force_index_compression, surface_position_index_compression, surface_position_data_compression, periodic_BC_data_compression, force_fabric_tensor_data_compression, kinetic_energy_data_compression = bertil_data_gatherer(
            simulation_directory + '_compression')
    elif simulation_directory.startswith("c:"):
        force_data_compression, surface_force_index_compression, surface_position_index_compression, surface_position_data_compression, periodic_BC_data_compression, force_fabric_tensor_data_compression, kinetic_energy_data_compression = local_data_gatherer(
            simulation_directory + '_compression')
    else:
        print("Error with simulation directory")

    if simulation_directory.startswith("/scratch"):
        force_data_tension, surface_force_index_tension, surface_position_index_tension, surface_position_data_tension, periodic_BC_data_tension, force_fabric_tensor_data_tension, kinetic_energy_data_tension = bertil_data_gatherer(
            simulation_directory + '_tension')
    elif simulation_directory.startswith("c:"):
        force_data_tension, surface_force_index_tension, surface_position_index_tension, surface_position_data_tension, periodic_BC_data_tension, force_fabric_tensor_data_tension, kinetic_energy_data_tension = local_data_gatherer(
            simulation_directory + '_tension')
    else:
        print("Error with simulation directory")
    time_tension, linear_strain_tension, sxx_tension, syy_tension, szz_tension = stress_and_linear_strain_finder(periodic_BC_data_tension, force_fabric_tensor_data_tension,surface_position_data_tension)
    time_compression, linear_strain_compression, sxx_compression, syy_compression, szz_compression = stress_and_linear_strain_finder(periodic_BC_data_compression, force_fabric_tensor_data_compression,surface_position_data_compression)

    # ==FIG 1 STRESS STRAIN IN COMPRESSION==============================================================================
    fig_compression,ax_compression = plt.subplots()
    ax_compression.set_ylabel('Stress [MPa]')
    ax_compression.set_xlabel('Strain [%]')
    lns_compression_xx = ax_compression.plot(linear_strain_compression[:] * 100, sxx_compression[:] / 1e6,
                                             label=r'$\sigma_{xx}$')
    lns_compression_yy = ax_compression.plot(linear_strain_compression[:] * 100, syy_compression[:] / 1e6,
                                             label=r'$\sigma_{yy}$')
    lns_compression_zz = ax_compression.plot(linear_strain_compression[:] * 100, szz_compression[:] / 1e6,
                                             label=r'$\sigma_{zz}$')
    ax_compression.set_title('Macroscopic stress in compression')
    plt.legend(loc='best')

    fname = fig_dir + 'compression_stress_strain'
    plt.savefig(fname)

    # ==FIG 2 STRESS STRAIN IN TENSION==================================================================================
    fig_tension,ax_tension = plt.subplots()
    ax_tension.set_ylabel('Stress [MPa]')
    ax_tension.set_xlabel('Strain [%]')
    lns_tension_xx = ax_tension.plot(linear_strain_tension[:] * 100, sxx_tension[:] / 1e6, label=r'$\sigma_{xx}$')
    lns__tension_yy = ax_tension.plot(linear_strain_tension[:] * 100, syy_tension[:] / 1e6, label=r'$\sigma_{yy}$')
    lns_tension_zz = ax_tension.plot(linear_strain_tension[:] * 100, szz_tension[:] / 1e6, label=r'$\sigma_{zz}$')
    ax_tension.set_title('Macroscopic stress in tension')
    plt.legend(loc='best')

    fname = fig_dir + 'tension_stress_strain'
    plt.savefig(fname)
    # ==FIG 3 STRESS AND STRAIN TO TIME FOR TENSION=====================================================================
    fig_tension_time, ax_tension_time = plt.subplots()
    ax_tension_time.set_ylabel("Stress [MPa]")
    ax_tension_time.set_xlabel("time [s]")
    lns_tension_xx_time = ax_tension_time.plot(time_tension, sxx_tension[:] / 1e6, label=r'$\sigma_{xx}$')
    lns_tension_yy_time = ax_tension_time.plot(time_tension, syy_tension[:] / 1e6,
                                                       label=r'$\sigma_{yy}$')
    lns_tension_zz_time = ax_tension_time.plot(time_tension, szz_tension[:] / 1e6,
                                                       label=r'$\sigma_{zz}$')
    ax_tension_time_2 = ax_tension_time.twinx()
    ##    ax2.set_xlabel("Calendering surface position [m]")
    #    ax2.set_ylabel("Calendering surface position [m]")
    lns_tension_time_2 = ax_tension_time_2.plot(time_tension, linear_strain_tension * 100, 'c', label=r'$\varepsilon_{xx}$')
    ax_tension_time_2.set_ylabel("Strain [%]")
    # lns3 = ax2.plot(time[max_index], linear_strain[max_index],'b',label='Position',linestyle="None",marker='x',markeredgewidth = 3)
    # lns4 = ax2.plot(time[min_index], linear_strain[min_index], 'b', label='Position', linestyle="None", marker='x',markeredgewidth = 3)

    lns_tension = lns_tension_xx_time + lns_tension_yy_time+ lns_tension_zz_time + lns_tension_time_2
    labs_tension = [l.get_label() for l in lns_tension]
    ax_tension_time.legend(lns_tension, labs_tension, loc=0)
    ax_tension_time.set_title('Macroscopic stress and strain in tension')

    fname = fig_dir + 'tension_stress_strain_time'
    plt.savefig(fname)

    # ==FIG 4 STRESS AND STRAIN TO TIME FOR COMPRESSION=================================================================
    fig_compression_time, ax_compression_time = plt.subplots()
    ax_compression_time.set_ylabel("Stress [MPa]")
    ax_compression_time.set_xlabel("time [s]")
    lns_compression_xx_time = ax_compression_time.plot(time_compression, sxx_compression[:] / 1e6, label=r'$\sigma_{xx}$')
    lns_compression_yy_time = ax_compression_time.plot(time_compression, syy_compression[:] / 1e6,
                                                       label=r'$\sigma_{yy}$')
    lns_compression_zz_time = ax_compression_time.plot(time_compression, szz_compression[:] / 1e6, label=r'$\sigma_{zz}$')
    ax_compression_time_2 = ax_compression_time.twinx()
    lns_compression_time_2 = ax_compression_time_2.plot(time_compression, linear_strain_compression * 100, 'c', label=r'$\varepsilon_{xx}$')
    ax_compression_time_2.set_ylabel("Strain [%]")
    lns_compression = lns_compression_xx_time + lns_compression_yy_time + lns_compression_zz_time + lns_compression_time_2
    labs_compression = [l.get_label() for l in lns_compression]
    ax_compression_time.legend(lns_compression, labs_compression, loc=0)
    ax_compression_time.set_title('Macroscopic stress and strain in compression')
    # ax_compression_time.set_ylim(ymin=-10,ymax=0)
    fig_compression_time.tight_layout()

    fname = fig_dir + 'compression_stress_strain_time'
    plt.savefig(fname)


    # ==FIG 5 KINETIC ENERGY IN TENSION=================================================================================
    fig_KE_tension,ax_KE_tension = plt.subplots()
    ax_KE_tension.set_ylabel('Kinetic energy [J]')
    ax_KE_tension.set_xlabel('Time [s]')
    lns_KE_tension = ax_KE_tension.plot(kinetic_energy_data_tension[: , -1].astype(float),
                                        kinetic_energy_data_tension[:,0].astype(float), label=r'KE')
    ax_KE_tension.set_title('Kinetic energy in tension')
    plt.legend(loc='best')

    fname = fig_dir + 'tension_KE'
    plt.savefig(fname)
    # ==FIG 6 KINETIC ENERGY IN COMPRESSION=============================================================================
    fig_KE_compression,ax_KE_compression = plt.subplots()
    ax_KE_compression.set_ylabel('Kinetic energy [J]')
    ax_KE_compression.set_xlabel('Time [s]')
    lns_KE_compression = ax_KE_compression.plot(kinetic_energy_data_compression[: , -1].astype(float), kinetic_energy_data_compression[:,0].astype(float), label=r'KE')
    ax_KE_compression.set_title('Kinetic energy in compression')
    plt.legend(loc='best')

    fname = fig_dir + 'compression_KE'
    plt.savefig(fname)

    # ==FIG 7 STIFFNESS AT STRAIN POINTS================================================================================
    if stiffness_at_points_flag == 1:
        strain_points_tension, stiffness_values_tension = stiffness_finder(sxx_tension, linear_strain_tension)
        strain_points_compression, stiffness_values_compression = stiffness_finder(sxx_compression,
                                                                                   linear_strain_compression)

        strain_points_compression = np.flip(strain_points_compression)
        stiffness_values_compression = np.flip(stiffness_values_compression)

        strain_points_total = np.concatenate((strain_points_compression, strain_points_tension))
        stiffness_values_total = np.concatenate((stiffness_values_compression, stiffness_values_tension))

        fig_stiff, ax_stiff = plt.subplots()
        ax_stiff.set_ylabel('Stiffness [GPa]')
        ax_stiff.set_xlabel('Strain [%]')
        ax_stiff.set_title('Stiffness of electrode layer')
        lns_stiff = ax_stiff.plot(strain_points_total * 100, stiffness_values_total * 1E-9, linestyle='dashed',
                                  marker='x',
                                  markerfacecolor='blue', markersize=12, markeredgewidth=3, linewidth=3)

        lns_stiff_exp_eps_dot_01 = ax_stiff.plot(exp_strain_points, Modulus_eps_dot_01,
                                                 marker='x', markersize=12, markeredgewidth=3,
                                                 linewidth=0, label='0.1 mm/min')

        lns_stiff_exp_eps_dot_05 = ax_stiff.plot(exp_strain_points, Modulus_eps_dot_05,
                                                 marker='x', markersize=12, markeredgewidth=3,
                                                 linewidth=0, label='0.5 mm/min')
        lns_stiff_exp_eps_dot_10 = ax_stiff.plot(exp_strain_points, Modulus_eps_dot_10,
                                                 marker='x', markersize=12, markeredgewidth=3,
                                                 linewidth=0, label='1.0 mm/min')
        lns_stiff_exp_eps_dot_100 = ax_stiff.plot(exp_strain_points, Modulus_eps_dot_100,
                                                  marker='x', markersize=12, markeredgewidth=3,
                                                  linewidth=0, label='10 mm/min')
        lns_stiff_exp_eps_dot_300 = ax_stiff.plot(exp_strain_points, Modulus_eps_dot_300,
                                                  marker='x', markersize=12, markeredgewidth=3,
                                                  linewidth=0, label='30 mm/min')

        plt.legend(loc='best')

        ax_stiff.set_ylim(ymin=0)
        fname = fig_dir + 'stiffness_points'
        plt.savefig(fname)
    # ==SHOW PLOT=======================================================================================================
    plt.show()
