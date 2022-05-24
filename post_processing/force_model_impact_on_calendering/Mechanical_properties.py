from calendering_pressure import local_data_gatherer,bertil_data_gatherer


import numpy as np

import matplotlib.pyplot as plt
import matplotlib
matplotlib.style.use('classic')




def main():

    simulation_directory = 'c:/Users/Axel/Documents/DEM/results/electrode_mechanical_loading/SN00_test'
    if simulation_directory.startswith("/scratch"):force_data, surface_force_index, surface_position_index, surface_position_data, periodic_BC_data, force_fabric_tensor_data = bertil_data_gatherer(simulation_directory)
    elif simulation_directory.startswith("c:"):force_data, surface_force_index, surface_position_index, surface_position_data, periodic_BC_data, force_fabric_tensor_data = local_data_gatherer(simulation_directory)
    else: print("Error with simulation directory")

    time = surface_position_data[:,-1]
    x_side_leght = periodic_BC_data[:, 2]-periodic_BC_data[:, 1]
    print(x_side_leght)
    x_side_leght_0 =  x_side_leght[0]
    y_side_length = periodic_BC_data[:, 4]-periodic_BC_data[:, 3]
    t0 = 1.11 # Chould be chaned later, defined in another way?=========================================================
    vol = x_side_leght*y_side_length*t0
    print(vol)

    sxx = force_fabric_tensor_data[:,1]/vol
    syy = force_fabric_tensor_data[:,5]/vol
    szz = force_fabric_tensor_data[:,9]/vol

    linear_strain = (x_side_leght[:]-x_side_leght_0)/x_side_leght_0

    plt.figure(0)
    plt.plot(linear_strain, sxx / 1e6, 'r', lw=2, label=r'$\sigma_{xx}$')
    plt.plot(linear_strain, syy / 1e6, 'g', lw=2, label=r'$\sigma_{yy}$')
    plt.plot(linear_strain, szz / 1e6, 'b', lw=2, label=r'$\sigma_{zz}$')
    plt.xlabel('Strain [-]')
    plt.ylabel('Stress [MPa]')
    plt.legend(loc='best')
    plt.show()

if __name__ == '__main__':
    main()