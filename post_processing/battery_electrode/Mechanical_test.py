
import numpy as np

import matplotlib.pyplot as plt
import matplotlib

matplotlib.style.use('classic')


def dimensions_box(data_directory):
    with open(data_directory + '/periodic_bc.dou', 'r') as periodic_bc:
        first_line = periodic_bc.readlines()[0]
        first_line = first_line.split(', ')
    id_idx = [i for i in range(len(first_line))]
    wall_data = np.genfromtxt(data_directory + '/periodic_bc.dou', delimiter=', ')
    data = np.zeros((wall_data.shape[0], 1))
    data = wall_data[:,  id_idx[0]+2]
    #print(data)
    return data


def position_zz(data_directory):
    with open(data_directory + '/surface_positions.dou', 'r') as periodic_bc:
        first_line = periodic_bc.readlines()[0]
        first_line = first_line.split(', ')
    id_idx = [i for i in range(len(first_line)) if first_line[i].startswith('ID')]
    surface_types = [first_line[idx+1][5:] for idx in id_idx]
    if surface_types.count('ID=') == 0 and surface_types.count('PointSurface') == 2:
        id_idx.sort(key=lambda x: first_line[x+1])
        wall_data = np.genfromtxt(data_directory + '/surface_positions.dou', delimiter=', ')
        time = np.zeros((wall_data.shape[0], 1))
        time = wall_data[:, id_idx[0]+32]

    return time



def pressures_box(data_directory):
    with open(data_directory + '/force_fabric_tensor.dou', 'r') as force_fabric_tensor:
        first_line = force_fabric_tensor.readlines()[0]
        first_line = first_line.split(', ')
    id_idx = [i for i in range(len(first_line))]
    force_data = np.genfromtxt(data_directory + '/force_fabric_tensor.dou', delimiter=', ')
    force = np.zeros((force_data.shape[0], 1))
    force = force_data[:,  id_idx[0]+1]
    return force

def pressures_box_yy(data_directory):
    with open(data_directory + '/force_fabric_tensor.dou', 'r') as force_fabric_tensor:
        first_line = force_fabric_tensor.readlines()[0]
        first_line = first_line.split(', ')
    id_idx = [i for i in range(len(first_line))]
    force_data = np.genfromtxt(data_directory + '/force_fabric_tensor.dou', delimiter=', ')
    force = np.zeros((force_data.shape[0], 1))
    force = force_data[:,  id_idx[0]+4]
    return force

def time_box(data_directory):
    with open(data_directory + '/force_fabric_tensor.dou', 'r') as force_fabric_tensor:
        first_line = force_fabric_tensor.readlines()[0]
        first_line = first_line.split(', ')
    id_idx = [i for i in range(len(first_line))]
    force_data = np.genfromtxt(data_directory + '/force_fabric_tensor.dou', delimiter=', ')
    time = np.zeros((force_data.shape[0], 1))
    time = force_data[:,  id_idx[0]]
    return time


if __name__ == '__main__':
    simulation_directory = '../../results/viscoelastic/4000-porosity_relaxation/'
    box_width = 0.726136 * 2
    box_height = 1.60915
    surface_height = 0.899473
    E = 2e9
    strain = (box_width - dimensions_box(simulation_directory)[5344:19081])/box_width
    Stress = pressures_box(simulation_directory)[5344:19081]/(box_width * surface_height * box_width *2)
    Stress_y = pressures_box_yy(simulation_directory)[5344:19081]/(box_width * box_height * box_width *2)
    stress = pressures_box(simulation_directory)[5344:19081]/(box_width * box_height *
                                                               dimensions_box(simulation_directory)[5344:19081] *2)

    #inkompresibelt på binder, isotropiskt material
    # inelastic strain
    epsilon_zz = -(position_zz(simulation_directory)[5344:19081]-surface_height)/surface_height
    print(position_zz(simulation_directory)[5344:19081])
    print(epsilon_zz)
    total_stress = Stress+Stress_y
    nu = -(E*epsilon_zz)/ total_stress
    print(nu)
    time = time_box(simulation_directory)[5344:19081]
    plt.plot(time, Stress)
    plt.xlabel("time[s]")
    plt.ylabel("Stress [Pa]")
    plt.show()
    plt.plot(strain, stress)
    plt.xlabel("Strain")
    plt.ylabel("Stress [Pa]")
    plt.show()
    plt.plot(time, Stress_y)
    plt.xlabel("Time")
    plt.ylabel("Stress_y[Pa]")

    plt.show()
    plt.plot(time, nu)
    plt.xlabel("Time")
    plt.ylabel("nu")

    plt.show()


