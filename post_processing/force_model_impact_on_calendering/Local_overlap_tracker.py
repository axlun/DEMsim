import numpy as np
import matplotlib.pyplot as plt
import os
import re
import pandas as pd
#from Local_particle_tracker import particle_tracker_local
plt.style.use('axel_style')


if __name__ == '__main__':
    # simulation_directory = 'c:/Users/Axel/Documents/DEM/Bertil_results/article_2/final_runs/particle_contact_model/SN_1_2/electrode_calendering_el_pl_binder_el_pl_particle/'
    # simulation_directory = "C:/Users/Axel/Documents/DEM/results/particle_tests/swelling/swelling_periodic_compaction/"
    # simulation_directory = "C:/Users/Axel/Documents/DEM/results/periodic_bc_tests/restart/"
    # simulation_directory = "C:/Users/Axel/Documents/DEM/results/periodic_bc_tests/no_restart/"
    # simulation_directory = 'c:/Users/Axel/Documents/DEM/results/electrode_natural_periodic_packing_el_pl_binder_el_pl_particle/SN_0/'
    simulation_directory = 'c:/Users/Axel/Documents/DEM/results/periodic_bc_tests/no_restart/'

    contact_files = os.listdir(simulation_directory + '/contacts')

    time = []
    contact_time_and_file_name_dict = {}
    for c_file in contact_files:
        contact_file = c_file
        time_stamp = re.split(r'\Acontacts_',contact_file)
        time_stamp = re.split(r'.dou\Z', time_stamp[1])

        time.append(float(time_stamp[0]))
        contact_time_and_file_name_dict[time_stamp[0]] = contact_file

    time.sort()
    time_vec = []
    overlap_data_vec = np.zeros((0, 0))
    particle_data_vec = []
    step = 1
    n_overlaps = 20
    overlap_matrix = np.zeros(shape=(int(len(time)/step)-1, 2 * n_overlaps))
    counter = 0
    for i in range(1, len(time), step):
        key = str(time[i])
        if time[i].is_integer():
            key = str(int(time[i]))
        file_to_open = simulation_directory+'contacts/' + contact_time_and_file_name_dict[key]
        data = pd.read_csv(file_to_open,header=None).to_numpy()
        sorted_overlap_data = np.sort(data[:,5])
        overlap_vector = np.zeros(2*n_overlaps)
        overlap_vector[:n_overlaps] = sorted_overlap_data[:n_overlaps]
        overlap_vector[n_overlaps:] = sorted_overlap_data[-n_overlaps:]
        # print(overlap_vector)
        # print(sorted_overlap_data)
        overlap_matrix[i-1,:] = overlap_vector
        # if i == 0:
        #     print(len(time),data[:,0].size)
        #     overlap_data_vec = np.zeros(shape=(data[:,0].size,int(len(time)/step)))
        #     overlap_data_vec[:,counter] = data[:,8]
        #     counter += 1
        #     particle_data_vec = data[np.where(data[:, 0] == p1)]
        # else:
        #     overlap_data_vec[:,counter] = data[:,8]
        #     counter += 1
            # particle_data_vec = np.vstack([particle_data_vec,data[np.where(data[:, 0] == p1)]])
        time_vec.append(float(key))
    time_vec = np.array(time_vec, dtype=float)


    figure_ke,ax_ke = plt.subplots()
    lns_ke = ax_ke.plot(time_vec, overlap_matrix)
    ax_ke.set_ylabel('Overlap [m]')
    ax_ke.set_xlabel('time [s]')

    plt.show()
