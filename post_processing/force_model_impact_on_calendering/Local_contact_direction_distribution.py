from Bertil_functions.Bertil_functions import *
from Bertil_calendering_pressure import local_data_gatherer, bertil_data_gatherer

import numpy as np

import matplotlib.pyplot as plt
import sys
import os
import re
import math
import matplotlib

matplotlib.style.use('axel_style')

def contact_direction_distribution_local(sim_dir):

    contact_files = os.listdir(sim_dir + '/contacts')

    # Gather all times for results files and sort them
    time = []
    contact_time_and_file_name_dict = {}
    for i in range(0,len(contact_files)):
        contact_file = contact_files[i]
        time_stamp = re.split(r'\Acontacts_',contact_file)
        time_stamp = re.split(r'.dou\Z', time_stamp[1])

        time.append(float(time_stamp[0]))
        contact_time_and_file_name_dict[time_stamp[0]] = contact_file

    time.sort()

    # Iterate trough all contact files and if particle or binder contact is activated, gather their z-normal component
    time_vec = []
    binder_contact_vec = []
    particle_contact_vec = []
    for i in range(0,len(time)):
        binder_contact = []
        particle_contact = []
        key = str(time[i])
        if time[i].is_integer():
            key = str(int(time[i]))
        ## Read file here
        file_to_open = sim_dir+'contacts/'+contact_time_and_file_name_dict[key]
        with open(file_to_open) as opened_contact_file:
            lines = opened_contact_file.readlines()
        # print(lines)

        for j in range(0, len(lines)):
            line_data = lines[j].split(', ')
            if float(line_data[6]) != 0 and float(line_data[8]) == 0.0:
                binder_contact.append(line_data[4])
            if float(line_data[6]) != 0 and float(line_data[8]) != 0:
                particle_contact.append(line_data[4])


        binder_contact_vec.append(binder_contact)
        particle_contact_vec.append(particle_contact)
        time_vec.append(float(key))

    time_vec = np.array(time_vec, dtype=float)
    p_length = max(map(len, particle_contact_vec))
    b_length = max(map(len, binder_contact_vec))

    particle_contact_mat = np.array([xi + [None]*(p_length-len(xi)) for xi in particle_contact_vec], dtype=float)
    binder_contact_mat = np.array([xi + [None]*(b_length-len(xi)) for xi in binder_contact_vec], dtype=float)


    return time_vec,particle_contact_mat,binder_contact_mat

if __name__ == '__main__':

    # simulation_directory = "C:/Users/Axel/Documents/DEM/results/particle_tests/swelling/swelling_periodic_compaction/no_restart/"
    simulation_directory = "C:/Users/Axel/Documents/DEM/Bertil_results/article_2/final_runs_2/SN_101/1/electrode_calendering_hertz/"
    # simulation_directory = "C:/Users/Axel/Desktop/contact_distribution_test/"
    time_vec,particle_contact_mat, binder_contact_mat = contact_direction_distribution_local(simulation_directory)

    for i in range(particle_contact_mat.shape[0]):
        for j in range(particle_contact_mat.shape[1]):
            ang = math.acos(particle_contact_mat[i,j])
            if ang > math.pi/2: ang = math.pi - ang
            particle_contact_mat[i,j] = ang

    n_bins = 10
    bin_size = math.pi/(2*n_bins)
    bin_edges = np.array([bin_size * i for i in range(n_bins+1)], dtype=float)
    bins = np.zeros((time_vec.size, n_bins))
    total_contacts = np.zeros((time_vec.size))
    bin_area = np.zeros(n_bins)
    for i in range(n_bins):
        bin_area[i] = 2 * math.pi * (math.cos(i*bin_size) - math.cos((i+1)*bin_size))
    for i in range(particle_contact_mat.shape[0]):
        for j in range(particle_contact_mat.shape[1]):
            if math.isnan(particle_contact_mat[i,j]) or particle_contact_mat[i, j] == 0:
                continue
            bin_no = int(particle_contact_mat[i,j] / bin_size)
            if bin_no == n_bins: bin_no -= 1
            bins[i, bin_no] += 1
            total_contacts[i] += 1
    particle_contacts_for_times = bins.sum(axis=1)
    counter = 0
    for i in range(0, len(time_vec), 20):
        if counter > 25: break
        counter += 1
        plt.figure()
        plt.stairs(bins[i,:]/(bin_area[:] * total_contacts[i]), bin_edges, fill=True)
        plt.xlabel('Angle to z-axis')
        plt.ylabel('Number of contacts')
        plt.title('Time t = ' + str(time_vec[i]))
        plt.ylim((0,1.2))
    print('show plots')
    plt.show()