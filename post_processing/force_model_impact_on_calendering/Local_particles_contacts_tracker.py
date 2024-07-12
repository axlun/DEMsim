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


def particles_contacts_tracker_local(sim_dir,p1):

    contact_files = os.listdir(sim_dir+'/contacts')
    time = []
    contact_time_and_file_name_dict = {}
    for i in range(0,len(contact_files)):
        contact_file = contact_files[i]
        time_stamp = re.split(r'\Acontacts_',contact_file)
        time_stamp = re.split(r'.dou\Z', time_stamp[1])

        time.append(float(time_stamp[0]))
        contact_time_and_file_name_dict[time_stamp[0]] = contact_file
    time.sort()
    time_vec = []
    contact_data_time_dict = {}
    contact_df_vec = []
    no_contacts = 0
    particles_in_contact = []
    for i in range(0,len(time)):
        contact_data_vec = []
        key = str(time[i])
        if time[i].is_integer():
            key = str(int(time[i]))
        ## Read file here
        file_to_open = sim_dir+'contacts/'+contact_time_and_file_name_dict[key]

        # print(file_to_open)
        data = pd.read_csv(file_to_open,header=None).to_numpy()
        df = pd.read_csv(file_to_open,header=None)
        # df = df.drop(df[df[0] != p1 | df[1] != p1].index)
        df = df.drop(df[(df[0] != p1) & (df[1] != p1)].index)
        if df.shape[0] > no_contacts: no_contacts = df.shape[0]
        for particle in df[0]:
            if particle not in particles_in_contact: particles_in_contact.append(particle)
        for particle in df[1]:
            if particle not in particles_in_contact: particles_in_contact.append(particle)
        contact_df_vec.append(df.to_numpy())
        # print(df)
        # if i == 1:
        #     contact_data_vec =np.zeros(len(data[0,:]))
        #     time_vec.append(float(key))
        #     continue
        # if len(data[np.where((data[:, 0] == p1))]) != 0:
        #     contact_data_vec = np.vstack([contact_data_vec,data[np.where((data[:, 0] == p1) )]])
        # if len(data[np.where((data[:, 1] == p1))]) != 0:
        #     contact_data_vec = np.vstack([contact_data_vec, data[np.where( (data[:, 1] == p1))]])
        # else:
        #     contact_data_vec = np.vstack([contact_data_vec,np.zeros(len(data[0,:]))])
        time_vec.append(float(key))
    particles_in_contact.remove(p1)
    particles_in_contact.sort()
    # print(particles_in_contact)
    contact_data_array = np.zeros((len(particles_in_contact), 20, len(time)))
    # print(contact_data_array.shape)
    for i in range(contact_data_array.shape[2]):
        contact_data_array[:, 0, i] = particles_in_contact
    # print(contact_data_array)

    for i in range(len(contact_df_vec)): # for all time steps in contact dataframe
        for j in range(len(contact_data_array[:, 0, 0])): # for all particles in contact array
            for k in range(len(contact_df_vec[i][:, 0])): # for all particles in data frame time step
                if contact_df_vec[i][k, 0] == contact_data_array[j, 0, i] or contact_df_vec[i][k,1] == contact_data_array[j,0,i]:
                    contact_data_array[j,1:,i] = contact_df_vec[i][k,2:]
    # print(contact_df_vec[0])
    # print(contact_data_array[:,:,0])
    time_vec = np.array(time_vec, dtype=float)
    return time_vec,contact_data_array

if __name__ == '__main__':
    # simulation_directory = 'c:/Users/Axel/Documents/DEM/results/swelling_electrode/SN_3/electrode_swelling/'
    # simulation_directory = 'c:/Users/Axel/Documents/DEM/results/swelling_electrode/SN_3/swelling_electrode_calendering/'
    # particle = 428

    # simulation_directory = 'c:/Users/Axel/Documents/DEM/results/swelling_electrode/SN_4/swelling_periodic_packing/'
    # particle = 569
    # simulation_directory = 'c:/Users/Axel/Documents/DEM/results/swelling_electrode/SN_6/electrode_swelling/'
    # particle = 447
    # particle = 32
    # simulation_directory = 'c:/Users/Axel/Documents/DEM/results/swelling_electrode/SN_5/swelling_electrode_calendering/'
    # particle = 241
    simulation_directory = 'c:/Users/Axel/Documents/DEM/results/swelling_electrode/SN_5/swelling_electrode_calendering/'
    particle = 427

    time_vec, contact_data_array = particles_contacts_tracker_local(simulation_directory,particle)

    fig_contact_force, ax_contact_force = plt.subplots()
    lns_total_contact = ax_contact_force.plot(time_vec,contact_data_array[:, 5, :].T) #, label=contact_data_array[:,0,0]) # , 'g', label=r'Total')
    # lns_binder_contacts = ax_contact_force.plot(time_vec, contact_data_vec[:,7], 'r', label=r'Binder')
    # lnd_particle_contacts = ax_contact_force.plot(time_vec,contact_data_vec[:,8],'b', label=r'Particle')
    ax_contact_force.set_ylabel("Force [N]")
    ax_contact_force.set_xlabel("time [s]")
    ax_contact_force.legend(loc='best')
    """
    fig_force_x_y_z, ax_force_x_y_z = plt.subplots()
    lns_normal_x = ax_force_x_y_z.plot(time_vec,contact_data_vec[:,6]*contact_data_vec[:,3], label='Normal x')
    lns_normal_y = ax_force_x_y_z.plot(time_vec,contact_data_vec[:,6]*contact_data_vec[:,4], label='Normal y')
    lns_normal_z = ax_force_x_y_z.plot(time_vec,contact_data_vec[:,6]*contact_data_vec[:,5], label='Normal z')
    lns_tangential_x = ax_force_x_y_z.plot(time_vec, contact_data_vec[:,10], label='Tangential x')
    lns_tangential_y = ax_force_x_y_z.plot(time_vec, contact_data_vec[:, 11], label='Tangential y')
    lns_tangential_z = ax_force_x_y_z.plot(time_vec, contact_data_vec[:, 12], label='Tangential z')

    ax_force_x_y_z.set_xlabel('Time [s]')
    ax_force_x_y_z.set_ylabel('Force [N]')
    ax_force_x_y_z.legend(loc='best')
    """

    fig_tangential_force, ax_tangential_force = plt.subplots()
    lns_tangential_force = ax_tangential_force.plot(time_vec, ((contact_data_array[: ,9, :]**2+
                                                                      contact_data_array[: ,10, :]**2+
                                                                      contact_data_array[: ,11, :]**2)**.5).T)
    ax_tangential_force.set_xlabel('Time [s]')
    ax_tangential_force.set_ylabel('Tangential force [N]')

    fig_overlap, ax_overlap = plt.subplots()
    lns_h = ax_overlap.plot(time_vec, contact_data_array[:, 4, :].T)
    ax_overlap.set_ylabel('Overlap [m]')
    ax_overlap.set_xlabel('Time [m]')

    fig_max_overlap, ax_max_overlap = plt.subplots()
    lns_h_max = ax_max_overlap.plot(time_vec,contact_data_array[:, 8, :].T)
    ax_max_overlap.set_ylabel("Max overlap [m]")
    ax_max_overlap.set_xlabel("time [s]")
    # ax_overlap.legend(loc='best')


    print('Show plot')
    plt.show()