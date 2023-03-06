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

    simulation_directory =  "C:/Users/Axel/Documents/DEM/results/contact_testing/elastic_plastic_binder_hertz_particle/wall_contact/meter_particle/"
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
        p1_data.append(data)
        contact_data.append(np.genfromtxt((simulation_directory+'contacts/'+particle_time_and_file_name_dict[key]).replace("particles","contacts"),delimiter=', '))

    force_fabric_tensor = np.genfromtxt(simulation_directory+'force_fabric_tensor.dou',delimiter=', ')
    print(contact_data)
    p1_data_mat = np.stack(p1_data,axis=0)
    contact_data_mat = np.stack(contact_data,axis=0)
    #print(force_fabric_tensor[1:,1]/(p1_data_mat[:,1]-p2_data_mat[:1]))

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

#=====================================OVERLAP TO CONTACT FORCE========================================================
    figure_overlap_force, ax_overlap_force = plt.subplots()
    lns_overlap_force = ax_overlap_force.plot(contact_data_mat[:,5],contact_data_mat[:,6],'r',linewidth=3,label=r'Contact force')
    lns_binder_force = ax_overlap_force.plot(contact_data_mat[:,5],contact_data_mat[:,7],'g--',linewidth=3,label=r'Binder force')
    lns_particle_force = ax_overlap_force.plot(contact_data_mat[:,5],contact_data_mat[:,8],'b--',linewidth=3,label=r'Particle force')
    ax_overlap_force.set_xlabel("Overlap [m]")
    ax_overlap_force.set_ylabel("Contact force [N]")
    ax_overlap_force.legend()

    plt.show()





    #
    #
    # matplotlib.style.use('classic')
    # plt.rc('font', serif='Computer Modern Roman')
    # plt.rcParams.update({'font.size': 15})
    # plt.rcParams['lines.linewidth'] = 2
    # force_data = np.genfromtxt('C:/Users/Axel/Documents/DEM/results/contact_testing/hertz_particle/contact_testing.dou',
    #                             delimiter=',')
    # simulation_dir = 'C:/Users/Axel/Documents/DEM/DEMsim/simulations/force_model_impact_on_electrode/contact_test.sim'
    #
    #
    # time_parameter = get_parameter(simulation_dir,'t')
    # R_eff = get_parameter(simulation_dir,'R')/2
    # Syp = get_parameter(simulation_dir,'particle_yield_stress_')
    # time = force_data[:,4]
    # ticks = time/time_parameter
    # plt.figure(1)
    # Ep = get_parameter(simulation_dir,'Ep')
    # nup = get_parameter(simulation_dir,'nup')
    # E0 = Ep/(1-nup**2)/2
    # R0 = 0.03/2
    # h = force_data[:, 0]
    # F = force_data[:,1]
    # h_norm = h/(R_eff)
    # F_norm = F/(R_eff**2 * Syp)
    # plt.plot(h, F)
    # plt.ylabel("Force [N]")
    # plt.xlabel("Overlap [m]")
    # plt.tight_layout()
    #
    # # print(ticks)
    # plt.figure(2)
    # plt.plot(ticks,F)
    # plt.ylabel("Force [N]")
    # plt.xlabel("Ticks [-]")
    # plt.tight_layout()
    #
    #
    # plt.figure(3)
    # plt.plot(time,F)
    # plt.ylabel("Force [N]")
    # plt.xlabel("Time [s]")
    # plt.tight_layout()
    #
    # plt.figure(4)
    # plt.plot(time,h)
    # plt.ylabel("Overlap [m]")
    # plt.xlabel("Time [s]")
    # plt.tight_layout()
    #
    #
    # plt.figure(5)
    # plt.plot(h_norm, F_norm)
    # plt.ylabel("Force/R²_eff [-]")
    # plt.xlabel("Overlap/R_eff  [-]")
    # plt.tight_layout()
    #
    # plt.show()
