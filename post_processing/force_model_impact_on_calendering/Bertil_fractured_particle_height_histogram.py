from Bertil_functions.Bertil_functions import *

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from io import StringIO
import pandas as pd
plt.style.use('axel_style')


if __name__ == '__main__':
    n_bins = 25
    n_bins = 12
    simulation_directory = '/scratch/users/axlun/DEMsim/results/article_4/particle_fracture/2/electrode_calendering_fracture_SF_1E6'
    time_step = 37.3639

    if (one_file_reader(simulation_directory+'/particles/particles_'+str(time_step)+'.dou') == OSError):
        print("Particle file not existing on Bertil")
        exit()

    particle_string = ''.join(one_file_reader(simulation_directory + '/particles/particles_' + str(time_step) + '.dou'))
    particle_file = StringIO(particle_string)
    particle_df = pd.read_csv(particle_file,
                              names=['Id','x','y','z','Rot_x','Rot_y','Rot_z','R','E_k_p','Mat-id','F_x', 'F_y','F_z'])
    particle_id_height_dict = {}
    print(particle_df['Id'].values)
    print(particle_df['Id'].size)
    for i in range(particle_df['Id'].size):
        particle_id_height_dict[particle_df['Id'].values[i]] = particle_df['z'].values[i]*100
    particle_heights = particle_df['z'].values*100

    simulation_time, number_of_fractures, fracture_array = fractured_particle_gatherer(simulation_directory)
    for row in fracture_array:
        if row[-1] == time_step:
            fracture_array = row[:-1]
    fractured_particle_height = np.zeros(np.size(fracture_array))

    counter = 0
    for id in fracture_array:
        fractured_particle_height[counter] = particle_id_height_dict[id]
        counter += 1

    hist_height, bin_edges_height = np.histogram(particle_heights, bins=n_bins)
    hist_height_fractured, bin_edges_height_fractured = np.histogram(fractured_particle_height, bins=bin_edges_height)

    figure_height, ax_height = plt.subplots()
    ax_height.stairs(hist_height,bin_edges_height,fill=True)
    # ax_height.set_title('Number of particles at height in electrode')
    ax_height.set_xlabel('Particle height')
    ax_height.set_ylabel('Number of particles')
    figure_height.tight_layout()

    figure_height_fractured, ax_height_fractured= plt.subplots()
    ax_height_fractured.stairs(hist_height_fractured,bin_edges_height_fractured,fill=True)
    # ax_height_fractured.set_title('Number of fractured particles at height in electrode')
    ax_height_fractured.set_xlabel('Particle height')
    ax_height_fractured.set_ylabel('Number of fractured particles')
    figure_height_fractured.tight_layout()

    hist_height_fractured_normalised = hist_height_fractured/hist_height

    figure_normalised_height_fractured, ax_normalised_height_fractured= plt.subplots()
    ax_normalised_height_fractured.stairs(hist_height_fractured_normalised,bin_edges_height,fill=True)
    # ax_normalised_height_fractured.set_title('Normalised number of fractured particles at height in electrode')
    ax_normalised_height_fractured.set_xlabel('Particle height')
    ax_normalised_height_fractured.set_ylabel('Normalised fractured particles to particle height')
    figure_normalised_height_fractured.tight_layout()

    plt.show()