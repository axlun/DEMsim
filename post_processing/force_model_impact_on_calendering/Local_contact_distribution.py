from Bertil_functions.Bertil_functions import *
from Bertil_calendering_pressure import local_data_gatherer, bertil_data_gatherer

import numpy as np

import matplotlib.pyplot as plt
import sys
import os
import re
import matplotlib

matplotlib.style.use('axel_style')

def contact_counter_local(sim_dir):

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
    binder_contact_vec = []
    particle_contact_vec = []
    binder_particle_contact_vec = []
    for i in range(0,len(time)):
        binder_contact = 0
        binder_particle_contact = 0
        particle_contact = 0
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
                binder_contact += 1
             if float(line_data[6]) != 0 and float(line_data[8]) != 0:
                particle_contact += 1
             if float(line_data[6]) != 0 and float(line_data[7]) != 0 and float(line_data[8]) != 0 and float(line_data[6]) == float(line_data[7]) + float(line_data[8]):
                 binder_particle_contact += 1


        binder_contact_vec.append(binder_contact)
        particle_contact_vec.append(particle_contact)
        binder_particle_contact_vec.append(particle_contact)
        time_vec.append(float(key))


    time_vec = np.array(time_vec, dtype=float)
    particle_contact_vec = np.array(particle_contact_vec, dtype=float)
    binder_contact_vec = np.array(binder_contact_vec, dtype=float)
    return time_vec,particle_contact_vec,binder_contact_vec, binder_particle_contact_vec

def kinetic_energy_local(sim_dir):

    kinetic_energy_data = np.genfromtxt(simulation_directory + '/kinetic_energy.dou', delimiter=',')
    time_vec = kinetic_energy_data[:,-1]
    ek_vec = kinetic_energy_data[:,0]
    #file_to_open = sim_dir + 'kinetic_energy.dou'
    # with open(file_to_open) as opened_contact_file:
    #     lines = opened_contact_file.readlines()
    # time_vec = []
    # ek_vec = []
    # for x in lines:
    #     time_vec.append(float(x.split(', ')[-1][:-1]))
    #     ek_vec.append(float(x.split(', ')[0]))




    return time_vec,ek_vec

if __name__ == '__main__':

    # simulation_directory = 'c:/Users/Axel/Documents/DEM/results/swelling_electrode/SN_5/swelling_periodic_packing/'
    simulation_directory = 'c:/Users/Axel/Documents/DEM/results/swelling_electrode/SN_5/swelling_electrode_calendering/'
    # simulation_directory = 'c:/Users/Axel/Documents/DEM/results/swelling_electrode/SN_5/swelling_electrode_mechanical_loading_compression/'
    # simulation_directory = 'c:/Users/Axel/Documents/DEM/results/swelling_electrode/SN_5/swelling_electrode_mechanical_loading_ss_0.9'
    # simulation_directory = 'c:/Users/Axel/Documents/DEM/results/swelling_electrode/SN_5/swelling_electrode_mechanical_loading_ss_0.95'

    time_vec, particle_contact_vec, binder_contact_vec, binder_particle_contact_vec = contact_counter_local(simulation_directory)
    time_vec_ek, ek_vec = kinetic_energy_local(simulation_directory)

    fig_contact_distribution, ax_contact_distribution = plt.subplots()
    lns_binder_contacts = ax_contact_distribution.plot(time_vec, binder_contact_vec, 'r', label=r'Binder contacts')
    lnd_particle_contacts = ax_contact_distribution.plot(time_vec,particle_contact_vec,'b', label=r'Particle contacts')
    ax_contact_distribution.set_ylabel("Number of contacts [-]")
    ax_contact_distribution.set_xlabel("time [s]")
    ax_contact_distribution.legend(loc='best')

    fig_kinetic_energy, ax_kinetic_energy = plt.subplots()
    lns_kinetic_energy = ax_kinetic_energy.plot(time_vec_ek, ek_vec, 'r', label=r'Kinetic energy')
    ax_kinetic_energy.set_ylabel("Kinetic energy [J]")
    ax_kinetic_energy.set_xlabel("time [s]")
    ax_kinetic_energy.legend(loc='best')
    plt.show()

