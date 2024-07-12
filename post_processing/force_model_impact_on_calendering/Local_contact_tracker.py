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

def contact_tracker_local(sim_dir,p1,p2):

    contact_files = os.listdir(sim_dir+'/contacts')
    particle_files = os.listdir(sim_dir+'/particles')
    time = []
    contact_time_and_file_name_dict = {}
    particle_time_and_file_name_dict = {}
    for i in range(0,len(contact_files)):
        contact_file = contact_files[i]
        particle_file = particle_files[i]
        time_stamp = re.split(r'\Acontacts_',contact_file)
        time_stamp = re.split(r'.dou\Z', time_stamp[1])

        time.append(float(time_stamp[0]))
        contact_time_and_file_name_dict[time_stamp[0]] = contact_file
        particle_time_and_file_name_dict[time_stamp[0]] = particle_file

    time.sort()
    time_vec = []
    contact_data_vec = []
    particle_overlap = []
    for i in range(0,len(time)):
        key = str(time[i])
        if time[i].is_integer():
            key = str(int(time[i]))

        #particle calc
        file_to_open = sim_dir+'particles/'+particle_time_and_file_name_dict[key]
        df = pd.read_csv(file_to_open,header=None)
        # df = df.drop(df[(df[0] != p1) | (df[0] != p2)].index).to_numpy()
        df = df.drop(df[~((df[0] == p1) | (df[0] == p2))].index).to_numpy()
        delta = ((df[0,1]-df[1,1])**2 + (df[0,2]-df[1,2])**2 + (df[0,3]-df[1,3])**2)**0.5
        r2 = df[0,7]+df[1,7]
        particle_overlap.append(r2 - delta)
        # Read file here
        file_to_open = sim_dir+'contacts/'+contact_time_and_file_name_dict[key]

        # print(file_to_open)
        data = pd.read_csv(file_to_open,header=None).to_numpy()
        if i == 0:
            if len(data[np.where((data[:, 0] == p1) * (data[:, 1] == p2))]) != 0:
                contact_data_vec = data[np.where((data[:, 0] == p1) * (data[:, 1] == p2))]
            elif len(data[np.where((data[:, 0] == p2) * (data[:, 1] == p1))]) != 0:
                contact_data_vec = data[np.where((data[:, 0] == p2) * (data[:, 1] == p1))]
            else:
                contact_data_vec = np.zeros(len(data[0, :]))
            # contact_data_vec =np.zeros(len(data[0,:]))
            time_vec.append(float(key))
            continue

        if len(data[np.where((data[:, 0] == p1) * (data[:, 1] == p2))]) != 0:
            contact_data_vec = np.vstack([contact_data_vec,data[np.where((data[:, 0] == p1) * (data[:, 1] == p2))]])
        elif len(data[np.where((data[:, 0] == p2) * (data[:, 1] == p1))]) != 0:
            contact_data_vec = np.vstack([contact_data_vec, data[np.where((data[:, 0] == p2) * (data[:, 1] == p1))]])
        else:
            contact_data_vec = np.vstack([contact_data_vec,np.zeros(len(data[0,:]))])
        time_vec.append(float(key))
    time_vec = np.array(time_vec, dtype=float)
    return time_vec,contact_data_vec, particle_overlap

if __name__ == '__main__':

    simulation_directory = 'c:/Users/Axel/Documents/DEM/results/swelling_electrode/SN_5/swelling_electrode_calendering/'
    particle_1 = 433
    particle_2 = 427

    time_vec, contact_data_vec, particle_overlap = contact_tracker_local(simulation_directory,particle_1,particle_2)

    fig_contact_force, ax_contact_force = plt.subplots()
    lns_total_contact = ax_contact_force.plot(time_vec,contact_data_vec[:,6], 'g', label=r'Total')
    lns_binder_contacts = ax_contact_force.plot(time_vec, contact_data_vec[:,7], 'r', label=r'Binder')
    lnd_particle_contacts = ax_contact_force.plot(time_vec,contact_data_vec[:,8],'b', label=r'Particle')
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

    fig_normal_tangential_force, ax_normal_tangential_force = plt.subplots()
    lns_normal_force = ax_normal_tangential_force.plot(time_vec, contact_data_vec[:,6], label='Normal')
    lns_tangential_force = ax_normal_tangential_force.plot(time_vec, (contact_data_vec[:,10]**2+
                                                                      contact_data_vec[:,11]**2+
                                                                      contact_data_vec[:,12]**2)**.5,
                                                           label='Tangential')
    ax_normal_tangential_force.set_xlabel('Time [s]')
    ax_normal_tangential_force.set_ylabel('Force [N]')
    ax_normal_tangential_force.legend(loc='best')

    fig_overlap, ax_overlap = plt.subplots()
    lns_h_ = ax_overlap.plot(time_vec,contact_data_vec[:,5], 'g', label=r'h_')
    lns_h_max = ax_overlap.plot(time_vec,contact_data_vec[:,9], 'r', label=r'h_max')
    lns_b_t_ = ax_overlap.plot(time_vec,-contact_data_vec[:,20], 'b', label=r'b_t')
    lns_calc_h = ax_overlap.plot(time_vec, particle_overlap, 'y', label=r'h_calc')
    # ax_calendering_surface_position = ax_calendering_surface_pressure.twinx()

    ax_binder_contact = ax_overlap.twinx()
    lns_binder_contact =  ax_binder_contact.plot(time_vec,contact_data_vec[:,19], 'c+', label=r'Binder contact')
    ax_binder_contact.set_ylabel("Binder contact [-]")
    ax_overlap.set_ylabel("Overlap [m]")
    ax_overlap.set_xlabel("time [s]")
    ax_overlap.legend(loc='best')

    lns = lns_h_ + lns_h_max + lns_b_t_ + lns_calc_h + lns_binder_contact
    labs = [l.get_label() for l in lns]
    ax_overlap.legend(lns, labs, loc='best')

    print('Show plot')
    plt.show()