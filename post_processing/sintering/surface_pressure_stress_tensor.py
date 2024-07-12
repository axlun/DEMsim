#
# Created by Axel  on 2024-05-28
#


import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
plt.style.use('axel_style')

# Intended for plotting compaction and swelling of particles in a cube with center in origo
if __name__=='__main__':
    # simulation_directory = "C:/Users/Axel/Documents/DEM/results/sintering/SN_1/"
    # simulation_directory = "C:/Users/Axel/Documents/DEM/results/sintering/rho_50/"
    # simulation_directory = "C:/Users/Axel/Documents/DEM/results/sintering/rho_60/"
    # simulation_directory = "C:/Users/Axel/Documents/DEM/results/sintering/rho_64/"
    # simulation_directory = "C:/Users/Axel/Documents/DEM/results/sintering/rho_64_5/"
    # simulation_directory = "C:/Users/Axel/Documents/DEM/results/sintering/rho_64_5/"
    # simulation_directory = "C:/Users/Axel/Documents/DEM/results/sintering/rho_64_6/"
    # simulation_directory = "C:/Users/Axel/Documents/DEM/results/sintering/rho_65/"
    # simulation_directory = "C:/Users/Axel/Documents/DEM/results/sintering/mixed/"
    # simulation_directory = "C:/Users/Axel/Documents/DEM/results/sintering/average/"
    # simulation_directory = "C:/Users/Axel/Documents/DEM/results/sintering/average_2/"
    # simulation_directory = "C:/Users/Axel/Documents/DEM/results/sintering/average_3/"
    # simulation_directory = "C:/Users/Axel/Documents/DEM/results/sintering/mixed_4/"
    # simulation_directory = "C:/Users/Axel/Documents/DEM/results/sintering/average_rho_64/"
    # simulation_directory = "C:/Users/Axel/Documents/DEM/results/sintering/mixed_rho_62/"
    simulation_directory = "C:/Users/Axel/Documents/DEM/results/sintering/mixed_rho_64_n_1e4/"
    # simulation_directory = "C:/Users/Axel/Documents/DEM/results/sintering/average_rho_64_n_1e4/"
    # simulation_directory = "C:/Users/Axel/Documents/DEM/results/sintering/average_rho_62/"


    surface_force_data = pd.read_csv(simulation_directory + '/surface_forces.dou', header=None).to_numpy()
    surface_position_data = pd.read_csv(simulation_directory + '/surface_positions.dou', header=None).to_numpy()
    force_fabric_tensor_data = pd.read_csv(simulation_directory + '/force_fabric_tensor.dou', header=None).to_numpy()

    simulation_time = surface_force_data[:,-1]
    cube_side = 2*surface_position_data[:,5]
    side_area = cube_side**2
    box_volume = cube_side**3

    stress_fabric_tensor = (force_fabric_tensor_data[:,1:].T/box_volume).T

    surface_forces = surface_force_data[:, 1:-1:5]
    surface_pressures = (surface_forces.T / side_area).T
    # print(surface_pressures)
    stress_fabric_tensor_directions = ['xx','xy','xz','yx','yy','yz','zx','zy','zz']
    fig_stress_fabric_tensor, ax_stress_fabric_tensor = plt.subplots()
    for i in range(9):
        ax_stress_fabric_tensor.plot(simulation_time, stress_fabric_tensor[:,i],
                                     label='$\sigma_{'+ stress_fabric_tensor_directions[i]+'}$')
    ax_stress_fabric_tensor.set_xlabel('Time [s]')
    ax_stress_fabric_tensor.set_ylabel('Stress [Pa]')
    ax_stress_fabric_tensor.legend(loc='best')


    fig_surface_pressure, ax_surface_pressure = plt.subplots()
    for i in range(6):
        ax_surface_pressure.plot(simulation_time, surface_pressures[:,i], label='surface ' + str(i))
    ax_surface_pressure.set_xlabel('Time [s]')
    ax_surface_pressure.set_ylabel('Pressure [Pa]')
    ax_surface_pressure.legend(loc='best')

    plt.show()




