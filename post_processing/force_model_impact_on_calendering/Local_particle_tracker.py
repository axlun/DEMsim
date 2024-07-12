from Bertil_functions.Bertil_functions import *
from Bertil_calendering_pressure import local_data_gatherer, bertil_data_gatherer

import numpy as np

import matplotlib.pyplot as plt
import sys
import os
import re
import matplotlib
import pandas as pd

matplotlib.style.use('axel_style')

def particle_tracker_local(sim_dir,p1):

    contact_files = os.listdir(sim_dir+'/contacts')
    particle_files = os.listdir(sim_dir+'/particles')
    time = []
    particle_time_and_file_name_dict = {}
    for i in range(0,len(particle_files)):
        particle_file = particle_files[i]
        time_stamp = re.split(r'\Aparticles_',particle_file)
        time_stamp = re.split(r'.dou\Z', time_stamp[1])

        time.append(float(time_stamp[0]))
        particle_time_and_file_name_dict[time_stamp[0]] = particle_file

    time.sort()
    time_vec = []
    particle_data_vec = []
    for i in range(0,len(time)):
        key = str(time[i])
        if time[i].is_integer():
            key = str(int(time[i]))
        ## Read file here
        file_to_open = sim_dir+'particles/'+particle_time_and_file_name_dict[key]

        data = pd.read_csv(file_to_open,header=None).to_numpy()

        if i == 0:
            particle_data_vec = data[np.where(data[:, 0] == p1)]
        else:
            particle_data_vec = np.vstack([particle_data_vec,data[np.where(data[:, 0] == p1)]])
        time_vec.append(float(key))
    time_vec = np.array(time_vec, dtype=float)
    return time_vec,particle_data_vec

if __name__ == '__main__':

    simulation_directory = 'c:/Users/Axel/Documents/DEM/results/swelling_electrode/SN_5/swelling_electrode_calendering/'
    particle = 427

    time_vec, particle_data_vec = particle_tracker_local(simulation_directory,particle)
    # print(time_vec)
    # print(particle_data_vec[:,0])
    # print(len(time_vec))
    # print(len(particle_data_vec))

    fig_particle_force_dir, ax_particle_force_dir = plt.subplots()
    lns_total_force_x = ax_particle_force_dir.plot(time_vec, particle_data_vec[:,10], label=r'x')
    lns_total_force_y = ax_particle_force_dir.plot(time_vec, particle_data_vec[:,11], label=r'y')
    lns_total_force_z = ax_particle_force_dir.plot(time_vec, particle_data_vec[:,12], label=r'z')
    ax_particle_force_dir.set_ylabel("Force [N]")
    ax_particle_force_dir.set_xlabel("time [s]")
    ax_particle_force_dir.legend(loc='best')

    fig_particle_force, ax_particle_force = plt.subplots()
    lns_total_force = ax_particle_force.plot(time_vec,(particle_data_vec[:,10]**2+particle_data_vec[:,11]**2+particle_data_vec[:,12]**2)**.5, 'g', label=r'Total')
    ax_particle_force.set_ylabel("Force [N]")
    ax_particle_force.set_xlabel("time [s]")
    ax_particle_force.legend(loc='best')

    fig_ke, ax_ke = plt.subplots()
    lns_h_ = ax_ke.plot(time_vec,particle_data_vec[:,8], 'r', label=r'KE')
    ax_ke.set_ylabel("Kinetic energy [J]")
    ax_ke.set_xlabel("time [s]")
    ax_ke.legend(loc='best')

    # =3D plot of particle position=====================================================================================
    # ax_position = plt.figure().add_subplots(projection='3d')
    ax_3D_pos = plt.figure().add_subplot(projection='3d')

    # ax_position = plt.figure().add_subplots(projection='3d')
    # fig_position,ax_position = plt.subplots(projection='3d')
    ax_3D_pos.plot(particle_data_vec[:,1],particle_data_vec[:,2],particle_data_vec[:,3],label='Particle position')
    ax_3D_pos.legend(loc='best')


    # =2D plot of particle position=====================================================================================
    fig_2D_pos, ax_2D_pos = plt.subplots()
    ax_2D_pos.scatter(particle_data_vec[:,1],particle_data_vec[:,2], c=range(0,len(particle_data_vec[:,1])))#colors#array form 0 to len(time_step))

    # =1D plot of particle positio==n===================================================================================
    fig_1D_pos, ax_1D_pos = plt.subplots()
    ax_1D_pos.plot(time_vec, particle_data_vec[:,1], label='X')
    ax_1D_pos.plot(time_vec, particle_data_vec[:,2], label='Y')
    ax_1D_pos.plot(time_vec, particle_data_vec[:,3], label='Z')
    ax_1D_pos.legend(loc='best')
    ax_1D_pos.set_xlabel('Time [s]')
    ax_1D_pos.set_ylabel('Position [m]')

    print('show plots')
    plt.show()