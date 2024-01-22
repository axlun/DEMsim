import numpy as np
import matplotlib.pyplot as plt
import os
import re
import pandas as pd

plt.style.use('axel_style')

if __name__ == '__main__':
    particle_file_directory = 'c:/Users/Axel/Documents/DEM/Bertil_results/article_2/final_runs/particle_contact_model/SN_1_2_2/electrode_mechanical_loading_el_pl_binder_el_pl_particle_compression/particles/particles_73.831.dou'
    # particle_file = os.listdir(particle_file_directory)
    particle_file_data = pd.read_csv(particle_file_directory,header=None).to_numpy()
    contact_counter = 0
    for i in range(0, (particle_file_data.shape[0]-1)):
        # print(particle_file_data[i,0])
        for j in range((i+1), particle_file_data.shape[0]):
            distance = ((particle_file_data[i, 1] - particle_file_data[j, 1]) ** 2 +
                        (particle_file_data[i, 2] - particle_file_data[j, 2]) ** 2 +
                        (particle_file_data[i, 3] - particle_file_data[j, 3]) ** 2) ** .5
            limit_distance = particle_file_data[i,7] + particle_file_data[j,7]
            if distance < limit_distance*.5:
                # print('contact')
                print('large overlap between particle: ' + str(int(particle_file_data[i,0])) + ' and ' +
                      str(int(particle_file_data[j,0])))
                print('particle distance is '+ str(distance))
                print('particle overlap is ' +str(limit_distance-distance))
                contact_counter += 1
    print(contact_counter)
    print('program end')