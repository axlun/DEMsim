import numpy as np
import matplotlib.pyplot as plt
import os
import re
import pandas as pd
#from Local_particle_tracker import particle_tracker_local
plt.style.use('axel_style')


if __name__ == '__main__':

    # simulation_directory = 'c:/Users/Axel/Documents/DEM/results/article_3/swelling/SN_2/swelling_electrode_calendering/'
    simulation_directory = 'c:/Users/Axel/Documents/DEM/results/article_3/swelling_material_scaling/SN_1/electrode_cycling/'

    particle_files = os.listdir(simulation_directory + '/particles')

    time = []
    particle_time_and_file_name_dict = {}
    n_particles = 0
    for p_file in particle_files:
        particle_file = p_file
        file_to_open = simulation_directory+'particles/' + p_file
        data = pd.read_csv(file_to_open,header=None).to_numpy()
        n_particles = data.shape[0]
        time_stamp = re.split(r'\Aparticles_',particle_file)
        time_stamp = re.split(r'.dou\Z', time_stamp[1])

        time.append(float(time_stamp[0]))
        particle_time_and_file_name_dict[time_stamp[0]] = particle_file

    time.sort()
    time_vec = []
    overlap_data_vec = np.zeros((0, 0))
    particle_data_vec = []
    step = 1
    # n_particles = 250
    height_matrix = np.zeros(shape=(int(len(time)/step)-1,  n_particles))
    n_in_avg = 100
    avg_height_vec = np.zeros(shape=(int(len(time)/step)-1))
    counter = 0
    for i in range(1, len(time), step):
        key = str(time[i])
        if time[i].is_integer():
            key = str(int(time[i]))
        file_to_open = simulation_directory+'particles/' + particle_time_and_file_name_dict[key]
        data = pd.read_csv(file_to_open,header=None).to_numpy()
        # sorted_particle_data = np.sort(data[:,3])
        particle_vector = np.zeros(n_particles)
        # particle_vector[:n_particles] = sorted_particle_data[:n_particles]
        # particle_vector[n_particles:] = sorted_particle_data[-n_particles:]  # Use this if lowest particles also
                                                                               # want to be captured
        particle_vector[:n_particles] = data[:,3] + data[:,7]
        height_matrix[i-1,:] = particle_vector
        # avg_height_vec[i-1] = np.sum(np.sort(particle_vector)[-n_in_avg:])/n_in_avg
        avg_height_vec[i-1] = np.sum(np.sort(data[:,3] + data[:,7])[-n_in_avg:])/n_in_avg
        time_vec.append(float(key))
    time_vec = np.array(time_vec, dtype=float)


    figure_ke,ax_ke = plt.subplots()
    lns_ke = ax_ke.plot(time_vec, height_matrix)
    lns_avg = ax_ke.plot(time_vec, avg_height_vec, '*')
    ax_ke.set_ylabel('Height [m]')
    ax_ke.set_xlabel('time [s]')

    figure_avg,ax_avg = plt.subplots()
    lns_avg = ax_avg.plot(time_vec, avg_height_vec, '*')
    ax_avg.set_ylabel('Height [m]')
    ax_avg.set_xlabel('time [s]')

    figure_avg_norm, ax_avg_norm = plt.subplots()
    lns_avg = ax_avg_norm.plot(time_vec, avg_height_vec/avg_height_vec[0], '*')
    ax_avg_norm.set_ylabel('Normalised height [m/m]')
    ax_avg_norm.set_xlabel('time [s]')
    plt.show()
