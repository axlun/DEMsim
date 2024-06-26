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

    simulation_directory = "C:/Users/Axel/Documents/DEM/results/contact_testing/elastic_plastic_binder_hertz_particle/tangential_force_wall_contact/New_tangential_force_relation1/"
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
        p1_data.append(data[0])

        contact_data.append(np.genfromtxt((simulation_directory+'contacts/'+particle_time_and_file_name_dict[key]).replace("particles","contacts"),delimiter=', '))

    force_fabric_tensor = np.genfromtxt(simulation_directory+'force_fabric_tensor.dou',delimiter=', ')
    kinetic_energy = np.genfromtxt(simulation_directory+'kinetic_energy.dou',delimiter=', ')

    p1_data_mat = np.stack(p1_data,axis=0)
    contact_data_mat = np.stack(contact_data,axis=0)
    #print(force_fabric_tensor[1:,1]/(p1_data_mat[:,1]-p2_data_mat[:1]))

    # =====================================PARTICLE 1 y POSITION========================================================
    figure_partcle_y_positions_time, ax_particle_y_positions_time = plt.subplots()
    lns_particle_1_y_pos_t = ax_particle_y_positions_time.plot(time, p1_data_mat[:,2],'r',linewidth=3,label=r'P_1')
    ax_particle_y_positions_time.set_xlabel('Time [s]')
    ax_particle_y_positions_time.set_ylabel('y position [m]')
    ax_particle_y_positions_time.legend()
# =====================================PARTICLE 1 X/Y POSITION========================================================
    figure_particle_positions, ax_particle_positions = plt.subplots()
    lns_particle_1 = ax_particle_positions.plot(p1_data_mat[:,1],p1_data_mat[:,2],'r*',linewidth=3,label=r'P_1')
    # lns_particle_2 = ax_particle_positions.plot(p2_data_mat[:,1],p2_data_mat[:,2],'b*',linewidth=3,label=r'P_2')
    ax_particle_positions.set_xlabel("x position [m]")
    ax_particle_positions.set_ylabel("y postition [m]")
    ax_particle_positions.legend()
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
    lns_overlap_force = ax_overlap_force.plot(contact_data_mat[:,5],contact_data_mat[:,9],'r',linewidth=3,label=r'Contact force')
    lns_binder_force = ax_overlap_force.plot(contact_data_mat[:,5],contact_data_mat[:,10],'g--',linewidth=3,label=r'Binder force')
    lns_particle_force = ax_overlap_force.plot(contact_data_mat[:,5],contact_data_mat[:,11],'b--',linewidth=3,label=r'Particle force')
    ax_overlap_force.set_xlabel("Overlap [m]")
    ax_overlap_force.set_ylabel("Contact force [N]")
    ax_overlap_force.legend()


#==KINETIC ENERGY=======================================================================================================
    figure_KE, ax_KE = plt.subplots()
    lns_KE = ax_KE.plot(kinetic_energy[:, -1], kinetic_energy[:, 0], 'r', linewidth=3,
                                              label=r'Kinetic Energy')
    ax_KE.set_xlabel("Time [S]")
    ax_KE.set_ylabel("Kinetic Energy [J]")
    ax_KE.legend()

#==CONTACT FORCE========================================================================================================

    fig_particle_forces, ax_particle_forces = plt.subplots()
    lns_particle_1_force = ax_particle_forces.plot(time,p1_data_mat[:,11],'r*-',linewidth=3,label=r'P_1')
    # lns_particle_2_force = ax_particle_forces.plot(time, p2_data_mat[:, 11], 'b', linewidth=3, label=r'P_2')
    #lns_force_fabric_tensor_xx = ax_particle_forces.plot(time,(force_fabric_tensor[:,1])/(p2_data_mat[:,1]-p1_data_mat[:,1]),'g--',linewidth=3,label=r'Force fabric XX')
    ax_particle_forces.set_xlabel("Simulation time [s]")
    ax_particle_forces.set_ylabel("Force in particle [N]")
    ax_particle_forces.legend()

    # ==CONTACT NORMAL ALL DIRECTIONS========================================================================================
    fig_all_contact_normals, ax_all_contact_normals = plt.subplots()
    lns_contact_1_normal_x = ax_all_contact_normals.plot(time, contact_data_mat[:, 2], 'r-', linewidth=3, label=r'n_x')
    lns_contact_1_normal_y = ax_all_contact_normals.plot(time, contact_data_mat[:, 3], 'g-', linewidth=3, label=r'n_y')
    lns_contact_1_normal_z = ax_all_contact_normals.plot(time, contact_data_mat[:, 4], 'b-', linewidth=3, label=r'n_z')
    ax_all_contact_normals.set_title("Normal of contact")
    ax_all_contact_normals.set_xlabel("Simulation time [s]")
    ax_all_contact_normals.set_ylabel("Normal direction [-]")
    ax_all_contact_normals.legend()
    fig_all_contact_normals.tight_layout()

    # ==CONTACT OVERLAP=====================================================================================================
    fig_contact_overlap, ax_contact_overlap = plt.subplots()
    lns_contact_overlap = ax_contact_overlap.plot(time, contact_data_mat[:, 5] + contact_data_mat[:, 20], 'r-',
                                                  linewidth=3, label=r'h_-b_t')
    ax_contact_overlap.set_title("Binder overlap")
    ax_contact_overlap.set_xlabel("Simulation time [s]")
    ax_contact_overlap.set_ylabel("Overlap [m]")
    ax_contact_overlap.legend()
    fig_contact_overlap.tight_layout()
# =====================================PARTICLE 1 X,Y,Z POSITION========================================================
    figure_partcle_positions_time, ax_particle_positions_time = plt.subplots()
    lns_particle_1_pos_x_t = ax_particle_positions_time.plot(time, p1_data_mat[:,1]+p1_data_mat[:,7]+contact_data_mat[0,-1],'r',linewidth=3,label=r'P_1_x')
    lns_particle_1_pos_y_t = ax_particle_positions_time.plot(time, p1_data_mat[:,2],'g',linewidth=3,label=r'P_1_y')
    lns_particle_1_pos_z_t = ax_particle_positions_time.plot(time, p1_data_mat[:,3],'b',linewidth=3,label=r'P_1_z')
    ax_particle_positions_time.set_title("Position of particle")
    ax_particle_positions_time.set_xlabel('Time [s]')
    ax_particle_positions_time.set_ylabel('Position [m]')
    ax_particle_positions_time.legend()
    figure_partcle_positions_time.tight_layout()


# ==PARTICLE FORCE ALL DIRECTIONS=======================================================================================
    fig_all_particle_forces, ax_all_particle_forces = plt.subplots()
    lns_particle_1_force_x = ax_all_particle_forces.plot(time, p1_data_mat[:, 10], 'r-', linewidth=3, label=r'F_x')
    lns_particle_1_force_y = ax_all_particle_forces.plot(time, p1_data_mat[:, 11], 'g-', linewidth=3, label=r'F_y')
    lns_particle_1_force_z = ax_all_particle_forces.plot(time, p1_data_mat[:, 12], 'b-', linewidth=3, label=r'F_z')
    lns_particle_1_force_res = ax_all_particle_forces.plot(time, (p1_data_mat[:, 10]**2 + p1_data_mat[:, 11]**2 + p1_data_mat[:, 12]**2)**(1/2), 'c--', linewidth=3, label=r'F_res')
    ax_all_particle_forces.set_title("Forces on particles")
    ax_all_particle_forces.set_xlabel("Simulation time [s]")
    ax_all_particle_forces.set_ylabel("Force on particle [N]")
    ax_all_particle_forces.legend()
    fig_all_particle_forces.tight_layout()

    # ==CONTACT FORCE ALL DIRECTIONS========================================================================================
    fig_all_contact_forces, ax_all_contact_forces = plt.subplots()
    lns_contact_1_force_n = ax_all_contact_forces.plot(time, contact_data_mat[:, 6], 'y-', linewidth=3, label=r'F_n')
    lns_contact_1_force_x = ax_all_contact_forces.plot(time, contact_data_mat[:, 10], 'r-', linewidth=3, label=r'F_Tx')
    lns_contact_1_force_y = ax_all_contact_forces.plot(time, contact_data_mat[:, 11], 'g-', linewidth=3, label=r'F_Ty')
    lns_contact_1_force_z = ax_all_contact_forces.plot(time, contact_data_mat[:, 12], 'b-', linewidth=3, label=r'F_Tz')
    lns_contact_1_force_res = ax_all_contact_forces.plot(time, (
                contact_data_mat[:, 10] ** 2 + contact_data_mat[:, 11] ** 2 + contact_data_mat[:, 12] ** 2) ** (1 / 2), 'c-',
                                                          linewidth=3, label=r'F_Tres')
    ax_all_contact_forces.set_title("Force in contact")
    ax_all_contact_forces.set_xlabel("Simulation time [s]")
    ax_all_contact_forces.set_ylabel("Force in particle [N]")
    ax_all_contact_forces.legend()
    fig_all_contact_forces.tight_layout()





    print('Showing Plot')
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
