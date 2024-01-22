import os

import numpy as np
import re
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.style


def get_parameter(parameter_dir, parameter):
    parameter_array = np.genfromtxt(parameter_dir, dtype='U30,f', delimiter='=')
    for index, item in np.ndenumerate(parameter_array):
        if parameter_array[index]['f0'] == parameter:
            return parameter_array[index]['f1']


if __name__ == '__main__':
    plt.style.use('axel_style')

    simulation_directory =  "C:/Users/Axel/Documents/DEM/results/contact_testing/elastic_plastic_binder_elastic_plastic_particle/wall_contact_drop/SN_2/"
    particle_files = os.listdir(simulation_directory+"particles/")
    time = []
    particle_time_and_file_name_dict = {}
    for i in range(0,len(particle_files)):
        contact_file = particle_files[i]
        time_stamp = re.split(r'\Aparticles_', contact_file)
        time_stamp = re.split(r'.dou\Z', time_stamp[1])
        time.append(float(time_stamp[0]))
        particle_time_and_file_name_dict[time_stamp[0]] = contact_file
    time.sort()
    p1_data = []
    overlap_data = []
    contact_data = []
    for i in range(0,len(time)):
        key = str(time[i])
        if time[i].is_integer():
            key = str(int(time[i]))
        # particle_file_to_open = simulation_directory+'particles/'+particle_time_and_file_name_dict[key]
        # with open(particle_file_to_open) as opened_contact_file:
        #     lines = opened_contact_file.readlines()
        # print(lines)
        # # print(simulation_directory+'particles/'+particle_time_and_file_name_dict[key])
        data = np.genfromtxt(simulation_directory+'particles/'+particle_time_and_file_name_dict[key],delimiter=', ')
        # print((data))
        p1_data.append(data)
        contact_data.append(
            np.genfromtxt((simulation_directory+'contacts/'+particle_time_and_file_name_dict[key]).replace(
                "particles","contacts"),delimiter=', '))
    p1_data_mat = np.stack(p1_data,axis=0)
    contact_data_mat = np.stack(contact_data, axis=0)
    # print(contact_data_mat)
    force_fabric_tensor = np.genfromtxt(simulation_directory+'force_fabric_tensor.dou',delimiter=', ')
    kinetic_energy = np.genfromtxt(simulation_directory+'kinetic_energy.dou',delimiter=', ')
    # print(p1_data_mat)
    kinetic_energy_particle_1 = p1_data_mat[:,7]
    potential_energy_particle_1 = 4/3 * 3.14* p1_data_mat[:,7] ** 3 * 4800 * 10 * (p1_data_mat[:,3] - p1_data_mat[:,7])
    #m*a*h

    # contact_data_mat = np.stack(contact_data,axis=0)

    fig_KE,ax_KE = plt.subplots()
    lns_ke = ax_KE.plot(kinetic_energy[:,-1],kinetic_energy[:,0],label='KE')
    lns_pe = ax_KE.plot(time,potential_energy_particle_1,label='PE')
    lns_totE = ax_KE.plot(kinetic_energy[:,-1][1:],potential_energy_particle_1+kinetic_energy[:,0][1:],label='Tot E')
    ax_KE.set_xlabel('Simulation time  [s]')
    ax_KE.set_ylabel('Kinetic energy [J]')
    ax_KE.legend(loc='best')


    fig_particle_forces, ax_particle_forces = plt.subplots()
    lns_particle_1_force = ax_particle_forces.plot(time,p1_data_mat[:,12],'r*-',linewidth=3,label=r'P_1')
    #lns_force_fabric_tensor_xx = ax_particle_forces.plot(time,(force_fabric_tensor[:,1])/(p2_data_mat[:,1]-p1_data_mat[:,1]),'g--',linewidth=3,label=r'Force fabric XX')
    ax_particle_forces.set_xlabel("Simulation time [s]")
    ax_particle_forces.set_ylabel("Z Force in particle [N]")
    ax_particle_forces.legend()


    figure_partcle_positions_time, ax_particle_positions_time = plt.subplots()
    lns_particle_1_pos_t = ax_particle_positions_time.plot(time, p1_data_mat[:,3],'r',linewidth=3,label=r'P_1')
    ax_particle_positions_time.set_xlabel('Time [s]')
    ax_particle_positions_time.set_ylabel('z position [m]')
    ax_particle_positions_time.legend()

    # figure_particle_positions, ax_particle_positions = plt.subplots()
    # lns_particle_1 = ax_particle_positions.plot(p1_data_mat[:,3],p1_data_mat[:,2],'r*',linewidth=3,label=r'P_1')
    # ax_particle_positions.set_xlabel("x position [m]")
    # ax_particle_positions.set_ylabel("y postition [m]")
    # ax_particle_positions.legend()
    #
    # lns5 = ax_sigma_zz_time2.plot(calendering_time, -szz*1e-6, 'r',linewidth=3,label=r'$\sigma_{zz}$')
    #
    # fig_full_calendering_sim, ax_full_calendering_sim = plt.subplots()
    # ax_full_calendering_sim.plot(1e2*calendering_surface_position[:], calendering_surface_pressure[:] * 1e-6)
    # ax_full_calendering_sim.set_ylabel("Calendering surface pressure [MPa]")
    # ax_full_calendering_sim.set_xlabel("Calendering surface position [µm]")
    # ax_full_calendering_sim.set_title('calendering surface pressure')

    # =====================================OVERLAP TO CONTACT FORCE=====================================================
    figure_overlap_force, ax_overlap_force = plt.subplots()
    lns_overlap_force = ax_overlap_force.plot(
        contact_data_mat[:,5],contact_data_mat[:,6],'r',linewidth=3,label=r'Contact force')
    lns_binder_force = ax_overlap_force.plot(
        contact_data_mat[:,5],contact_data_mat[:,7],'g--',linewidth=3,label=r'Binder force')
    lns_particle_force = ax_overlap_force.plot(
        contact_data_mat[:,5],contact_data_mat[:,8],'b--*',linewidth=3,label=r'Particle force')
    ax_overlap_force.set_xlabel("Overlap [m]")
    ax_overlap_force.set_ylabel("Contact force [N]")
    ax_overlap_force.legend()

    # =====================================OVERLAP TIME=================================================================
    figure_overlap_time, ax_overlap_time = plt.subplots()
    lns_overlap_time = ax_overlap_time.plot(time, contact_data_mat[:, 5], 'r', linewidth=3,
                                              label=r'overlap')
    # lns_binder_force = ax_overlap_force.plot(contact_data_mat[:, 5], contact_data_mat[:, 7], 'g--', linewidth=3,
    #                                         label=r'Binder force')
    # lns_particle_force = ax_overlap_force.plot(contact_data_mat[:, 5], contact_data_mat[:, 8], 'b--', linewidth=3,
    #                                           label=r'Particle force')
    ax_overlap_time.set_ylabel("Overlap [m]")
    ax_overlap_time.set_xlabel("Time [s]")
    ax_overlap_time.legend()

    plt.show()