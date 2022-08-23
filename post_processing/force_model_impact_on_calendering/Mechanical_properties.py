from calendering_pressure import local_data_gatherer, bertil_data_gatherer

import numpy as np

import matplotlib.pyplot as plt
import matplotlib

matplotlib.style.use('classic')


def stress_and_linear_strain_finder(periodic_BC_data, force_fabric_tensor_data):
    time = periodic_BC_data[:, 0]

    x_side_length = periodic_BC_data[:, 2] - periodic_BC_data[:, 1]
    # print(x_side_length)
    x_side_length_0 = x_side_length[0]
    y_side_length = periodic_BC_data[:, 4] - periodic_BC_data[:, 3]
    t0 = 1.11  # Chould be chaged later, defined in another way?=========================================================
    vol = x_side_length * y_side_length * t0
    # print(vol)

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
        ##==================STRAIN VÄRDEN LÄGGS INTE TILL I ARRAYN...===========
        stiffness_values = np.append(stiffness_values, (stress_vec[max_index[i]] - stress_vec[min_index[i]]) / (
                    strain_vec[max_index[i]] - strain_vec[min_index[i]]))
    if strain_vec[1] < 0:
        strain_points = strain_vec[min_index]
    if strain_vec[1] > 0:
        strain_points = strain_vec[max_index]

    return strain_points, stiffness_values  # returns strain in [-] and stiffness in Pa


if __name__ == '__main__':

    #    simulation_directory = 'c:/Users/Axel/Documents/DEM/results/electrode_mechanical_loading/SN00_test'
    simulation_directory = '/scratch/users/axlun/DEMsim/results/electrode_mechanical_loading/SN_hertz_5000p_btr_065_hal_85'
    if simulation_directory.startswith("/scratch"):
        force_data_compression, surface_force_index_compression, surface_position_index_compression, surface_position_data_compression, periodic_BC_data_compression, force_fabric_tensor_data_compression = bertil_data_gatherer(
            simulation_directory + '_compression')
    elif simulation_directory.startswith("c:"):
        force_data_compression, surface_force_index_compression, surface_position_index_compression, surface_position_data_compression, periodic_BC_data_compression, force_fabric_tensor_data_compression = local_data_gatherer(
            simulation_directory + '_compression')
    else:
        print("Error with simulation directory")

    if simulation_directory.startswith("/scratch"):
        force_data_tension, surface_force_index_tension, surface_position_index_tension, surface_position_data_tension, periodic_BC_data_tension, force_fabric_tensor_data_tension = bertil_data_gatherer(
            simulation_directory + '_tension')
    elif simulation_directory.startswith("c:"):
        force_data_tension, surface_force_index_tension, surface_position_index_tension, surface_position_data_tension, periodic_BC_data_tension, force_fabric_tensor_data_tension = local_data_gatherer(
            simulation_directory + '_tension')
    else:
        print("Error with simulation directory")

    time_tension, linear_strain_tension, sxx_tension, syy_tension, szz_tension = stress_and_linear_strain_finder(
        periodic_BC_data_tension, force_fabric_tensor_data_tension)
    strain_points_tension, stiffness_values_tension = stiffness_finder(sxx_tension, linear_strain_tension)

    time_compression, linear_strain_compression, sxx_compression, syy_compression, szz_compression = stress_and_linear_strain_finder(
        periodic_BC_data_compression, force_fabric_tensor_data_compression)
    strain_points_compression, stiffness_values_compression = stiffness_finder(sxx_compression,
                                                                               linear_strain_compression)



    fig_compression,ax_compression = plt.subplots()
    ax_compression.set_ylabel('Stress [MPa]')
    ax_compression.set_xlabel('Strain [%]')
    lns_compression_xx = ax_compression.plot(linear_strain_compression[:] * 100, sxx_compression[:] / 1e6, 'r', lw=2, label=r'$\sigma_{xx}$')
    lns_compression_yy = ax_compression.plot(linear_strain_compression[:] * 100, syy_compression[:] / 1e6, 'g', lw=2, label=r'$\sigma_{yy}$')
    lns_compression_zz = ax_compression.plot(linear_strain_compression[:] * 100, szz_compression[:] / 1e6, 'b', lw=2, label=r'$\sigma_{zz}$')
    ax_compression.set_title('calendering surface pressure in compression')
    plt.legend(loc='best')

    fig_tension,ax_tension = plt.subplots()
    ax_tension.set_ylabel('Stress [MPa]')
    ax_tension.set_xlabel('Strain [%]')
    lns_tension_xx = ax_tension.plot(linear_strain_tension[:] * 100, sxx_tension[:] / 1e6, 'r', lw=2, label=r'$\sigma_{xx}$')
    lns__tension_yy = ax_tension.plot(linear_strain_tension[:] * 100, syy_tension[:] / 1e6, 'g', lw=2, label=r'$\sigma_{yy}$')
    lns_tension_zz = ax_tension.plot(linear_strain_tension[:] * 100, szz_tension[:] / 1e6, 'b', lw=2, label=r'$\sigma_{zz}$')
    ax_tension.set_title('calendering surface pressure in tension')
    plt.legend(loc='best')




    fig_tension_time, ax_tension_time = plt.subplots()
    ax_tension_time.set_ylabel("Stress [MPa]")
    ax_tension_time.set_xlabel("time [s]")
    lns_tension_time = ax_tension_time.plot(time_tension, sxx_tension[:] / 1e6, 'r',lw=2, label='Pressure')
    ax_tension_time_2 = ax_tension_time.twinx()
    ##    ax2.set_xlabel("Calendering surface position [m]")
    #    ax2.set_ylabel("Calendering surface position [m]")
    lns_tension_time_2 = ax_tension_time_2.plot(time_tension, linear_strain_tension, 'b',lw=2, label='Position')
    ax_tension_time_2.set_ylabel("Strain [%]")
    # lns3 = ax2.plot(time[max_index], linear_strain[max_index],'b',label='Position',linestyle="None",marker='x',markeredgewidth = 3)
    # lns4 = ax2.plot(time[min_index], linear_strain[min_index], 'b', label='Position', linestyle="None", marker='x',markeredgewidth = 3)

    # added these three lines
    lns_tension = lns_tension_time + lns_tension_time_2
    labs_tension = [l.get_label() for l in lns_tension]
    ax_tension_time.legend(lns_tension, labs_tension, loc=0)
    ax_tension_time.set_title('calendering surface pressure')



    strain_points_compression = np.flip(strain_points_compression)
    stiffness_values_compression = np.flip(stiffness_values_compression)

    strain_points_total = np.concatenate((strain_points_compression, strain_points_tension))
    stiffness_values_total = np.concatenate((stiffness_values_compression, stiffness_values_tension))






    fig_stiff, ax_stiff = plt.subplots()
    ax_stiff.set_ylabel('Stiffness [GPa]')
    ax_stiff.set_xlabel('Strain [%]')
    lns_stiff = ax_stiff.plot(strain_points_total * 100, stiffness_values_total * 1E-9, color='red', linestyle='dashed',
                              marker='x',
                              markerfacecolor='blue', markersize=12, markeredgewidth=3, linewidth=3, label='bt = XXX')
    ax_stiff.set_ylim(ymin=0)
    plt.legend(loc='best')
    plt.show()
