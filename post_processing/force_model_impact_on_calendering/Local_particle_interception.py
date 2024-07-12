import numpy as np
import matplotlib.pyplot as plt
import os
import re
import pandas as pd

plt.style.use('axel_style')

if __name__ == '__main__':

    # ==SINTERING=======================================================================================================
    # simulation_directory = "C:/Users/Axel/Documents/DEM/results/sintering/SN_1/"

    # ==PERIODIC_PACKING================================================================================================
    simulation_directory = 'c:/Users/Axel/Documents/DEM/results/swelling_electrode/SN_5/swelling_periodic_packing/'

    particle_files = os.listdir(simulation_directory + '/particles')

    counter = 0
    for p_file in particle_files:
        # print(p_file)
        particle_file_data = pd.read_csv(simulation_directory + "/particles/" + p_file,header=None).to_numpy()
        mp_path = simulation_directory + "/mirror_particles/mirror_" + p_file
        if os.path.exists(mp_path) and os.path.getsize(mp_path) > 0:
            mirror_particle_file_data = pd.read_csv(mp_path, header=None).to_numpy()
            particle_file_data = np.concatenate((particle_file_data,mirror_particle_file_data), axis=0)
        time_stamp = re.split(r'\Aparticles_',p_file)
        time_stamp = re.split(r'.dou\Z', time_stamp[1])[0]
        counter += 1
        if counter == 10:
            print(time_stamp)
            counter = 0
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
                    print('particle distance is ' + str(distance))
                    print('particle overlap is ' + str(limit_distance-distance))
                    contact_counter += 1
        if contact_counter != 0:
            print(str(contact_counter) + " particles have large overlap at time t = " + str(time_stamp))
    print('program end')