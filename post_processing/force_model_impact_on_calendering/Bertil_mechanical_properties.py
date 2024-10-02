from force_model_impact_on_calendering.Bertil_calendering_pressure import local_data_gatherer, bertil_data_gatherer, contact_counter_bertil
from force_model_impact_on_calendering.Local_contact_distribution import contact_counter_local

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
    print('x and y side length is = ' + str(x_side_length[0]) + ', ' + str(y_side_length[0]))

    #t0 = 1.11  # Chould be chaged later, defined in another way?======================================================
    t0 = float(surface_position_data[0,32])-1
    print('layer height = ' + str(t0))
    vol = x_side_length * y_side_length * t0
    sxx = -force_fabric_tensor_data[:, 1] / vol
    syy = -force_fabric_tensor_data[:, 5] / vol
    szz = -force_fabric_tensor_data[:, 9] / vol

    tau_xy = -force_fabric_tensor_data[:, 2] / vol
    tau_xz = -force_fabric_tensor_data[:, 3] / vol

    tau_yx = -force_fabric_tensor_data[:, 4] / vol
    tau_yz = -force_fabric_tensor_data[:, 6] / vol

    tau_zx = -force_fabric_tensor_data[:, 7] / vol
    tau_zy = -force_fabric_tensor_data[:, 8] / vol


    linear_strain = (x_side_length[:] - x_side_length_0) / x_side_length_0

    return time, linear_strain, sxx, syy, szz, tau_xy, tau_xz, tau_yx, tau_yz, tau_zx, tau_zy


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


def mechanical_properties_plotting_func(simulation_directory, stiffness_at_points_flag=1, contact_flag=0):

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
    plt.style.use('axel_style')

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
        force_data_compression, surface_force_index_compression, surface_position_index_compression,\
            surface_position_data_compression, periodic_BC_data_compression, force_fabric_tensor_data_compression, \
            kinetic_energy_data_compression = bertil_data_gatherer(simulation_directory + '_compression')
        if contact_flag ==1:
            time_vec_compression, particle_contact_vec_compression, binder_contact_vec_compression, \
                binder_particle_contact_vec_compression = contact_counter_bertil(simulation_directory+ '_compression')
    elif simulation_directory.startswith("c:"):
        force_data_compression, surface_force_index_compression, surface_position_index_compression, \
        surface_position_data_compression, periodic_BC_data_compression, force_fabric_tensor_data_compression, \
        kinetic_energy_data_compression = local_data_gatherer(simulation_directory + '_compression')
        if contact_flag == 1:
            time_vec_compression, particle_contact_vec_compression, binder_contact_vec_compression, \
                binder_particle_contact_vec_compression = contact_counter_local(simulation_directory+ '_compression/')
    else:
        print("Error with simulation directory")

    if simulation_directory.startswith("/scratch"):
        force_data_tension, surface_force_index_tension, surface_position_index_tension, surface_position_data_tension, periodic_BC_data_tension, force_fabric_tensor_data_tension, kinetic_energy_data_tension = bertil_data_gatherer(
            simulation_directory + '_tension')
        if contact_flag == 1:
            time_vec_tension, particle_contact_vec_tension, binder_contact_vec_tension, binder_particle_contact_vec_tension = contact_counter_bertil(
                simulation_directory + '_tension')
    elif simulation_directory.startswith("c:"):
        force_data_tension, surface_force_index_tension, surface_position_index_tension, surface_position_data_tension, periodic_BC_data_tension, force_fabric_tensor_data_tension, kinetic_energy_data_tension = local_data_gatherer(
            simulation_directory + '_tension')
        if contact_flag == 1:
            time_vec_tension, particle_contact_vec_tension, binder_contact_vec_tension, binder_particle_contact_vec_tension = contact_counter_local(
                simulation_directory+ '_tension/')
    else:
        print("Error with simulation directory")
    time_tension, linear_strain_tension, sxx_tension, syy_tension, szz_tension, tau_xy_tension, tau_xz_tension, tau_yx_tension, tau_yz_tension, tau_zx_tension, tau_zy_tension = stress_and_linear_strain_finder(periodic_BC_data_tension, force_fabric_tensor_data_tension,surface_position_data_tension)
    time_compression, linear_strain_compression, sxx_compression, syy_compression, szz_compression, tau_xy_compression, tau_xz_compression, tau_yx_compression, tau_yz_compression, tau_zx_compression, tau_zy_compression = stress_and_linear_strain_finder(periodic_BC_data_compression, force_fabric_tensor_data_compression,surface_position_data_compression)

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
    lns_compression_tau_xy = ax_compression.plot(linear_strain_compression[:] * 100, tau_xy_compression[:] / 1e6, label=r'$\tau_{xy}$')
    lns_compression_tau_xz = ax_compression.plot(linear_strain_compression[:] * 100, tau_xz_compression[:] / 1e6, label=r'$\tau_{xz}$')
    lns_compression_tau_yx = ax_compression.plot(linear_strain_compression[:] * 100, tau_yx_compression[:] / 1e6, label=r'$\tau_{yx}$')
    lns_compression_tau_yz = ax_compression.plot(linear_strain_compression[:] * 100, tau_yz_compression[:] / 1e6, label=r'$\tau_{yz}$')
    lns_compression_tau_zx = ax_compression.plot(linear_strain_compression[:] * 100, tau_zx_compression[:] / 1e6, label=r'$\tau_{zx}$')
    lns_compression_tau_zy = ax_compression.plot(linear_strain_compression[:] * 100, tau_zy_compression[:] / 1e6, label=r'$\tau_{zy}$')


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

    lns_tension_tau_xy = ax_tension.plot(linear_strain_tension[:] * 100, tau_xy_tension[:] / 1e6, label=r'$\tau_{xy}$')
    lns_tension_tau_xz = ax_tension.plot(linear_strain_tension[:] * 100, tau_xz_tension[:] / 1e6, label=r'$\tau_{xz}$')
    lns_tension_tau_yx = ax_tension.plot(linear_strain_tension[:] * 100, tau_yx_tension[:] / 1e6, label=r'$\tau_{yx}$')
    lns_tension_tau_yz = ax_tension.plot(linear_strain_tension[:] * 100, tau_yz_tension[:] / 1e6, label=r'$\tau_{yz}$')
    lns_tension_tau_zx = ax_tension.plot(linear_strain_tension[:] * 100, tau_zx_tension[:] / 1e6, label=r'$\tau_{zx}$')
    lns_tension_tau_zy = ax_tension.plot(linear_strain_tension[:] * 100, tau_zy_tension[:] / 1e6, label=r'$\tau_{zy}$')


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
    # ==FIG 8 CONTACTS FOR STRAINS=======================================================================================================

    if (contact_flag == 1):
        rm_list_compression = []
        for i in range(len(time_compression)-1):
            if time_compression[i] == time_compression[i+1]:
                rm_list_compression.append(i)
        time_compression_short = np.delete(time_compression,rm_list_compression)
        linear_strain_compression_short = np.delete(linear_strain_compression,rm_list_compression)

        rm_list_tension = []
        for i in range(len(time_tension)-1):
            if time_tension[i] == time_tension[i+1]:
                rm_list_tension.append(i)
        time_tension_short = np.delete(time_tension,rm_list_tension)
        linear_strain_tension_short = np.delete(linear_strain_tension,rm_list_tension)

        fig_contacts_compression,ax_contacts_compression = plt.subplots()
        lns_contacts_compression = ax_contacts_compression.plot(linear_strain_compression_short[:300] * 100,
                                                                particle_contact_vec_compression[:300],'C0')
        lns_contacts_tension = ax_contacts_compression.plot(linear_strain_tension_short[:300] * 100,
                                                                particle_contact_vec_tension[:300],'C0')
        ax_contacts_compression.set_ylim(ymin=0)
        ax_contacts_compression.set_xlim(xmin=-2.2, xmax=2.2)
        ax_contacts_compression.set_ylabel('Particle contacts [-]')
        ax_contacts_compression.set_xlabel('Strain [%]')
        fname = fig_dir + 'contacts_to_strain'
        plt.savefig(fname)

if __name__ == '__main__':
    # simulation_directory = '/scratch/users/axlun/DEMsim/results/article_2/final_runs_2/SN_101/1/electrode_mechanical_loading_hertz'
    # simulation_directory = '/scratch/users/axlun/DEMsim/results/article_2/final_runs_2/SN_201_periodic_packing/3/electrode_mechanical_loading_el_pl_binder_el_pl_particle'
    # simulation_directory = 'c:/Users/Axel/Documents/DEM/results/swelling_electrode/SN_5/swelling_electrode_mechanical_loading'
    # simulation_directory = 'c:/Users/Axel/Documents/DEM/results/swelling_electrode/SN_5/swelling_electrode_mechanical_loading_ss_0.9'
    # simulation_directory = 'c:/Users/Axel/Documents/DEM/results/swelling_electrode/SN_5/swelling_electrode_mechanical_loading_ss_0.95'

    simulation_directory = '/scratch/users/axlun/DEMsim/results/article_3/charge_cycling/SN_1/swelling_electrode_mechanical_loading'
    # simulation_directory = '/scratch/users/axlun/DEMsim/results/article_3/charge_cycling/SN_1/swelling_electrode_mechanical_loading_material_scaling'
    # simulation_directory = '/scratch/users/axlun/DEMsim/results/article_3/charge_cycling/SN_1/swelling_electrode_mechanical_loading_material_scaling_05'
    # simulation_directory = '/scratch/users/axlun/DEMsim/results/article_3/charge_cycling/SN_1/swelling_electrode_mechanical_loading_ss_1.06266_material_scaling'
    # simulation_directory = '/scratch/users/axlun/DEMsim/results/article_3/charge_cycling/SN_1/swelling_electrode_mechanical_loading_ss_1.03228_material_scaling'
    # simulation_directory = '/scratch/users/axlun/DEMsim/results/article_3/charge_cycling/SN_1/swelling_electrode_mechanical_loading_ss_1.06266'
    # simulation_directory = '/scratch/users/axlun/DEMsim/results/article_3/charge_cycling/SN_1/swelling_electrode_mechanical_loading_cycle_1'
    # simulation_directory = '/scratch/users/axlun/DEMsim/results/article_3/charge_cycling/SN_1/swelling_electrode_mechanical_loading_cycle_3'
    # simulation_directory = '/scratch/users/axlun/DEMsim/results/article_3/charge_cycling/SN_1/swelling_electrode_mechanical_loading_cycle_10'

    # simulation_directory = '/scratch/users/axlun/DEMsim/results/article_3/charge_cycling/SN_1_reduced_cal/swelling_electrode_mechanical_loading_1135'
    # simulation_directory = '/scratch/users/axlun/DEMsim/results/article_3/charge_cycling/SN_1_reduced_cal/swelling_electrode_mechanical_loading_1115'


    # simulation_directory = '/scratch/users/axlun/DEMsim/results/article_3/charge_cycling/SN_2/swelling_electrode_mechanical_loading'
    # simulation_directory = '/scratch/users/axlun/DEMsim/results/article_3/charge_cycling/SN_2/swelling_electrode_mechanical_loading_material_scaling'
    # simulation_directory = '/scratch/users/axlun/DEMsim/results/article_3/charge_cycling/SN_2/swelling_electrode_mechanical_loading_ss_1.06266_material_scaling'
    # simulation_directory = '/scratch/users/axlun/DEMsim/results/article_3/charge_cycling/SN_2/swelling_electrode_mechanical_loading_ss_1.03228_material_scaling'
    # simulation_directory = '/scratch/users/axlun/DEMsim/results/article_3/charge_cycling/SN_2/swelling_electrode_mechanical_loading_ss_1.06266'
    # simulation_directory = '/scratch/users/axlun/DEMsim/results/article_3/charge_cycling/SN_2/swelling_electrode_mechanical_loading_cycle_1'
    # simulation_directory = '/scratch/users/axlun/DEMsim/results/article_3/charge_cycling/SN_2/swelling_electrode_mechanical_loading_cycle_3'
    simulation_directory = '/scratch/users/axlun/DEMsim/results/article_3/charge_cycling/SN_2/swelling_electrode_mechanical_loading_cycle_10'

    # simulation_directory = '/scratch/users/axlun/DEMsim/results/article_3/charge_cycling/SN_3/swelling_electrode_mechanical_loading'
    # simulation_directory = '/scratch/users/axlun/DEMsim/results/article_3/charge_cycling/SN_3/swelling_electrode_mechanical_loading_material_scaling'
    # simulation_directory = '/scratch/users/axlun/DEMsim/results/article_3/charge_cycling/SN_3/swelling_electrode_mechanical_loading_ss_1.06266_material_scaling'
    # simulation_directory = '/scratch/users/axlun/DEMsim/results/article_3/charge_cycling/SN_3/swelling_electrode_mechanical_loading_ss_1.03228_material_scaling'
    # simulation_directory = '/scratch/users/axlun/DEMsim/results/article_3/charge_cycling/SN_3/swelling_electrode_mechanical_loading_ss_1.06266'
    # simulation_directory = '/scratch/users/axlun/DEMsim/results/article_3/charge_cycling/SN_3/swelling_electrode_mechanical_loading_cycle_1'
    # simulation_directory = '/scratch/users/axlun/DEMsim/results/article_3/charge_cycling/SN_3/swelling_electrode_mechanical_loading_cycle_3'
    # simulation_directory = '/scratch/users/axlun/DEMsim/results/article_3/charge_cycling/SN_3/swelling_electrode_mechanical_loading_cycle_10'

    stiffness_at_points_flag = 1
    contact_flag = 0
    mechanical_properties_plotting_func(simulation_directory, stiffness_at_points_flag, contact_flag)
    # ==SHOW PLOT=======================================================================================================
    print("Plotting results")
    plt.show()
