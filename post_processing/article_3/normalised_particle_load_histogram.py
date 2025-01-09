import numpy as np
np.set_printoptions(threshold=np.inf)
import matplotlib
import matplotlib.pyplot as plt
from os.path import exists
import shutil
import os
from io import StringIO
import pandas as pd


def particle_load_normaliser(sim_dir, sim_time):
    contact_file = sim_dir + 'contacts/contacts_' + str(sim_time) + '.dou'
    particle_file = sim_dir + 'particles/particles_' + str(sim_time) + '.dou'
    particle_radius_dict = {}
    particle_fracture_dict = {}
    particle_df = pd.read_csv(particle_file, header=None)
    for index, row in particle_df.iterrows():
        particle_radius_dict[int(row[0])] = row[7]
        particle_fracture_dict[int(row[0])] = False
    contact_df = pd.read_csv(contact_file, header=None)
    normalised_particle_load_vec = np.zeros(contact_df.shape[0] * 2)
    particle_index_normalised_particle_load_matrix = np.zeros((contact_df.shape[0] * 2, 2))
    for index, row in contact_df.iterrows():
        if row[0] < 2 or row[1] < 2 or row[19] == 1:
            continue
        d1 = 2 * particle_radius_dict[int(row[0])]
        d2 = 2 * particle_radius_dict[int(row[1])]
        f = row[6]
        s1 = f / (d1 ** 2) * 1E-6
        s2 = f / (d2 ** 2) * 1E-6
        normalised_particle_load_vec[index*2] = s1
        normalised_particle_load_vec[index*2 + 1] = s2

        particle_index_normalised_particle_load_matrix[index*2, 0] = row[0]
        particle_index_normalised_particle_load_matrix[index*2, 1] = s1

        particle_index_normalised_particle_load_matrix[index*2+1, 0] = row[1]
        particle_index_normalised_particle_load_matrix[index*2+1, 1] = s2

        s_max = 100
        if s1 >= s_max:
            particle_fracture_dict[row[0]] = True
        if s2 >= s_max:
            particle_fracture_dict[row[1]] = True

    particle_fracture_array = np.array(list(particle_fracture_dict.items()))*2
    print(particle_fracture_array)
    plt.figure()
    plt.hist2d(particle_fracture_array[:,0], particle_fracture_array[:,1],bins=[5000,2])

    plt.figure()
    plt.hist(particle_fracture_array[:,1], bins=[0,1,2], density=True)

    normalised_particle_load_vec = normalised_particle_load_vec[normalised_particle_load_vec[:] != 0]
    particle_index_normalised_particle_load_matrix = \
        particle_index_normalised_particle_load_matrix[particle_index_normalised_particle_load_matrix[:,0] != 0]
    # plt.hist(normalised_particle_load_vec, bins=500)
    # print(particle_index_normalised_particle_load_matrix[:,0])
    # print(particle_index_normalised_particle_load_matrix[:,1])
    # plt.hist(particle_index_normalised_particle_load_matrix[:,1], bins=500)

    print(particle_fracture_dict)
    plt.figure()
    plt.hist2d(particle_index_normalised_particle_load_matrix[:,0],particle_index_normalised_particle_load_matrix[:,1],
               bins=[np.arange(0,5000,1),[0,150,1000]])
               # bins=[5000,15])
    plt.colorbar()

    plt.show()


if __name__ == '__main__':
    plt.style.use('axel_style')

    simulation_directory = 'c:/Users/Axel/Documents/DEM/Bertil_results/article_3/final_runs/1/swelling_electrode_calendering/'
    # time = 18.164
    # time = 23.564
    # time = 24.664
    time = 25.564
    # time = 28.564
    particle_load_normaliser(simulation_directory, time)
