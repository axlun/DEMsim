from Bertil_functions.Bertil_functions import *

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from io import StringIO
import pandas as pd
plt.style.use('axel_style')


if __name__ == '__main__':
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
    particle_id_size_dict = {}
    print(particle_df['Id'].values)
    print(particle_df['Id'].size)
    for i in range(particle_df['Id'].size):
        particle_id_size_dict[particle_df['Id'].values[i]] = particle_df['R'].values[i]*200
    particle_sizes = particle_df['R'].values*200

    simulation_time, number_of_fractures, fracture_array = fractured_particle_gatherer(simulation_directory)
    for row in fracture_array:
        if row[-1] == time_step:
            fracture_array = row[:-1]
    fractured_particle_sizes = np.zeros(np.size(fracture_array))

    counter = 0
    for id in fracture_array:
        fractured_particle_sizes[counter] = particle_id_size_dict[id]
        counter += 1

    hist_psd, bin_edges_psd = np.histogram(particle_sizes, bins=n_bins)
    hist_fractured, bin_edges_fractured = np.histogram(fractured_particle_sizes, bins=bin_edges_psd)

    figure_psd, ax_psd = plt.subplots()
    ax_psd.stairs(hist_psd,bin_edges_psd,fill=True)
    ax_psd.set_xlabel('Particle size [µm]')
    ax_psd.set_ylabel('Number of particles')
    figure_psd.tight_layout()

    figure_fractured, ax_fractured= plt.subplots()
    ax_fractured.stairs(hist_fractured,bin_edges_psd,fill=True)
    ax_fractured.set_xlabel('Particle size [µm]')
    ax_fractured.set_ylabel('Number of fractured particles')
    figure_fractured.tight_layout()

    hist_fractured_normalised = hist_fractured/hist_psd

    figure_normalised_fractured, ax_normalised_fractured= plt.subplots()
    ax_normalised_fractured.stairs(hist_fractured_normalised,bin_edges_psd,fill=True)
    ax_normalised_fractured.set_xlabel('Particle size [µm]')
    ax_normalised_fractured.set_ylabel('Normalised fractured particles to size')
    figure_normalised_fractured.tight_layout()

    plt.show()