import numpy as np
import matplotlib.pyplot as plt
import os
import re
import pandas as pd
#from Local_particle_tracker import particle_tracker_local
plt.style.use('axel_style')


if __name__ == '__main__':
    # simulation_directory = 'c:/Users/Axel/Documents/DEM/Bertil_results/article_2/final_runs/particle_contact_model/SN_1_2/electrode_calendering_el_pl_binder_el_pl_particle/'

    # simulation_directory = 'c:/Users/Axel/Documents/DEM/results/swelling_electrode/SN_5/swelling_electrode_calendering/'
    # simulation_directory = 'c:/Users/Axel/Documents/DEM/results/swelling_electrode/SN_5/electrode_swelling/'
    # simulation_directory = 'c:/Users/Axel/Documents/DEM/results/swelling_electrode/SN_6/electrode_swelling/'
    # simulation_directory = 'c:/Users/Axel/Documents/DEM/results/swelling_electrode/SN_6/swelling_electrode_mechanical_loading_tension/'
    # simulation_directory = 'c:/Users/Axel/Documents/DEM/results/swelling_electrode/SN_6/swelling_electrode_mechanical_loading_compression/'
    #simulation_directory = 'c:/Users/Axel/Documents/DEM/results/swelling_electrode/SN_6/swelling_electrode_mechanical_loading_ss_0.9_tension/'
    # simulation_directory = 'c:/Users/Axel/Documents/DEM/results/swelling_electrode/SN_6/swelling_electrode_mechanical_loading_ss_0.95_compression/'

    simulation_directory = "c:/Users/Axel/Documents/DEM/results/article_3/swelling_material_scaling/SN_1/electrode_cycling_1/"

    data = pd.read_csv(simulation_directory + "kinetic_energy.dou", header=None).to_numpy()
    total_time = data[:,-1]
    total_ke = data[:,2]

    figure_tot_ke,ax_tot_ke = plt.subplots()
    lns_ke = ax_tot_ke.plot(total_time, total_ke)
    ax_tot_ke.set_ylabel('Kinetic energy [J]')
    ax_tot_ke.set_xlabel('time [s]')

    particle_files = os.listdir(simulation_directory + '/particles')

    time = []
    particle_time_and_file_name_dict = {}
    for x in particle_files:
        particle_file = x
        time_stamp = re.split(r'\Aparticles_',particle_file)
        time_stamp = re.split(r'.dou\Z', time_stamp[1])

        time.append(float(time_stamp[0]))
        particle_time_and_file_name_dict[time_stamp[0]] = particle_file

    time.sort()
    time_vec = []
    ke_data_vec = np.zeros((0,0))
    radius_data_vec = np.zeros((0,0))
    particle_data_vec = []
    step = 1
    counter = 0
    for i in range(0, len(time), step):
        key = str(time[i])
        if time[i].is_integer():
            key = str(int(time[i]))
        ## Read file here
        file_to_open = simulation_directory+'particles/'+particle_time_and_file_name_dict[key]
        data = pd.read_csv(file_to_open,header=None).to_numpy()
        if i == 0:
            print(len(time),data[:,0].size)
            ke_data_vec = np.zeros(shape=(data[:,0].size,int(len(time)/step)))
            ke_data_vec[:,counter] = data[:,8]
            radius_data_vec = np.zeros(shape=(data[:,0].size,int(len(time)/step)))
            radius_data_vec[:,counter] = data[:,7]
            counter += 1
            # particle_data_vec = data[np.where(data[:, 0] == p1)]
        else:
            ke_data_vec[:,counter] = data[:,8]
            radius_data_vec[:,counter] = data[:,7]
            counter += 1
            # particle_data_vec = np.vstack([particle_data_vec,data[np.where(data[:, 0] == p1)]])
        time_vec.append(float(key))
    time_vec = np.array(time_vec, dtype=float)


    figure_ke,ax_ke = plt.subplots()
    lns_ke = ax_ke.plot(time_vec, ke_data_vec.T)
    ax_ke.set_ylabel('Kinetic energy [J]')
    ax_ke.set_xlabel('time [s]')

    figure_radius ,ax_radius = plt.subplots()
    lns_radius = ax_radius.plot(time_vec, radius_data_vec.T)
    ax_radius.set_ylabel('Particle radius [m]')
    ax_radius.set_xlabel('time [s]')

    plt.show()
