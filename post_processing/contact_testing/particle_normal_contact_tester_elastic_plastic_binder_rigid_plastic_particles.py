import os

import numpy as np
import re
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.style
plt.style.use('axel_style')


def get_parameter(parameter_dir, parameter):
    parameter_array = np.genfromtxt(parameter_dir, dtype='U30,f', delimiter='=')
    for index, item in np.ndenumerate(parameter_array):
        if parameter_array[index]['f0'] == parameter:
            return parameter_array[index]['f1']


if __name__ == '__main__':

    # simulation_directory =  "C:/Users/Axel/Documents/DEM/results/contact_testing/force_overlap/el_pl_particle/SN_1/gT"
    simulation_directory = "C:/Users/Axel/Documents/DEM/results/particle_tests/swelling/SN_1/"
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
    p2_data = []
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
        p2_data.append(data[1])
        contact_data.append(
            np.genfromtxt((simulation_directory+'contacts/'+particle_time_and_file_name_dict[key]).replace(
                "particles","contacts"),delimiter=', '))

    force_fabric_tensor = np.genfromtxt(simulation_directory+'force_fabric_tensor.dou',delimiter=', ')
    kinetic_energy = np.genfromtxt(simulation_directory+'kinetic_energy.dou',delimiter=', ')

    p2_data_mat = np.stack(p2_data,axis=0)
    p1_data_mat = np.stack(p1_data,axis=0)
    contact_data_mat = np.stack(contact_data,axis=0)
      #print(force_fabric_tensor[1:,1]/(p1_data_mat[:,1]-p2_data_mat[:1]))
    print(p1_data_mat[0,:])
    R_0 = 1/(1/p1_data_mat[:,7] + 1/p2_data_mat[:,7])
    print(R_0)

    # =====================================PARTICLE 1 y POSITION========================================================
    figure_partcle_y_positions_time, ax_particle_y_positions_time = plt.subplots()
    lns_particle_1_y_pos_t = ax_particle_y_positions_time.plot(time, p1_data_mat[:,2],'r',linewidth=3,label=r'P_1')
    ax_particle_y_positions_time.set_xlabel('Time [s]')
    ax_particle_y_positions_time.set_ylabel('y position [m]')
    ax_particle_y_positions_time.legend()

    # =====================================PARTICLE 1 X/Y POSITION======================================================
    figure_particle_positions, ax_particle_positions = plt.subplots()
    lns_particle_1 = ax_particle_positions.plot(p1_data_mat[:,1],p1_data_mat[:,2],'r*',linewidth=3,label=r'P_1')
    lns_particle_2 = ax_particle_positions.plot(p2_data_mat[:,1],p2_data_mat[:,2],'b*',linewidth=3,label=r'P_2')
    ax_particle_positions.set_xlabel("x position [m]")
    ax_particle_positions.set_ylabel("y position [m]")
    ax_particle_positions.legend()

    # =====================================OVERLAP TO CONTACT FORCE=====================================================
    figure_overlap_force, ax_overlap_force = plt.subplots()
    lns_overlap_force = ax_overlap_force.plot(contact_data_mat[:,5],contact_data_mat[:,6],'r',linewidth=3,label=r'Contact force')
    lns_binder_force = ax_overlap_force.plot(contact_data_mat[:,5],contact_data_mat[:,7],'g--',linewidth=3,label=r'Binder force')
    lns_particle_force = ax_overlap_force.plot(contact_data_mat[:,5],contact_data_mat[:,8],'b--',linewidth=3,label=r'Particle force')
    ax_overlap_force.set_xlabel("Overlap [m]")
    ax_overlap_force.set_ylabel("Contact force [N]")
    ax_overlap_force.legend()

    # =====================================NORMALISED OVERLAP TO CONTACT FORCE==========================================
    figure_norm_overlap_force, ax_norm_overlap_force = plt.subplots()
    lns_norm_overlap_force = ax_norm_overlap_force.plot(contact_data_mat[:,5]/R_0,contact_data_mat[:,6]/(R_0**2 * 4e9),
                                                        'r',linewidth=3,label=r'Contact force')
    lns_norm_binder_force = ax_norm_overlap_force.plot(contact_data_mat[:,5]/R_0,contact_data_mat[:,7]/(R_0**2 * 4e9),
                                                       'g--',linewidth=3,label=r'Binder force')
    lns_norm_particle_force = ax_norm_overlap_force.plot(contact_data_mat[:,5]/R_0,contact_data_mat[:,8]/(R_0**2 * 4e9),
                                                         'b--',linewidth=3,label=r'Particle force')
    ax_norm_overlap_force.set_xlabel("$h / R_0$")
    ax_norm_overlap_force.set_ylabel("$F / R_0^2 \sigma_Y$")
    ax_norm_overlap_force.legend()

    # ==KINETIC ENERGY==================================================================================================
    figure_KE, ax_KE = plt.subplots()
    lns_KE = ax_KE.plot(kinetic_energy[:, -1], kinetic_energy[:, 0], 'r', linewidth=3, label=r'Kinetic Energy')
    ax_KE.set_xlabel("Time [S]")
    ax_KE.set_ylabel("Kinetic Energy [J]")
    ax_KE.legend()

    # ==CONTACT FORCE========================================================================================================
    fig_particle_forces, ax_particle_forces = plt.subplots()
    lns_particle_1_force = ax_particle_forces.plot(time,p1_data_mat[:,10],'r*-',linewidth=3,label=r'P_1')
    lns_particle_2_force = ax_particle_forces.plot(time, p2_data_mat[:, 10], 'b', linewidth=3, label=r'P_2')
    # lns_force_fabric_tensor_xx = \
    #     ax_particle_forces.plot(time,(force_fabric_tensor[:,1])/(p2_data_mat[:,1]-p1_data_mat[:,1]),'g--',
    #     linewidth=3,label=r'Force fabric XX')
    ax_particle_forces.set_xlabel("Simulation time [s]")
    ax_particle_forces.set_ylabel("Force in x-dir in particle [N]")
    ax_particle_forces.legend()

    # ==CONTACT NORMAL ALL DIRECTIONS===================================================================================
    fig_all_contact_normals, ax_all_contact_normals = plt.subplots()
    lns_contact_1_normal_x = ax_all_contact_normals.plot(time, contact_data_mat[:, 2], 'r-', linewidth=3, label=r'n_x')
    lns_contact_1_normal_y = ax_all_contact_normals.plot(time, contact_data_mat[:, 3], 'g-', linewidth=3, label=r'n_y')
    lns_contact_1_normal_z = ax_all_contact_normals.plot(time, contact_data_mat[:, 4], 'b-', linewidth=3, label=r'n_z')
    ax_all_contact_normals.set_title("Normal of contact")
    ax_all_contact_normals.set_xlabel("Simulation time [s]")
    ax_all_contact_normals.set_ylabel("Normal direction [-]")
    ax_all_contact_normals.legend()
    fig_all_contact_normals.tight_layout()

    # ==CONTACT OVERLAP=================================================================================================
    fig_contact_overlap, ax_contact_overlap = plt.subplots()
    lns_contact_overlap = ax_contact_overlap.plot(time, contact_data_mat[:, 5] , 'r-', linewidth=3, label=r'h_')
    ax_contact_overlap.set_title("Particle overlap")
    ax_contact_overlap.set_xlabel("Simulation time [s]")
    ax_contact_overlap.set_ylabel("Overlap [m]")
    ax_contact_overlap.legend()
    fig_contact_overlap.tight_layout()

# =====================================PARTICLE 1 X,Y,Z POSITION========================================================
    figure_partcle_positions_time, ax_particle_positions_time = plt.subplots()
    lns_particle_1_pos_x_t = ax_particle_positions_time.plot(time, p1_data_mat[:,1], 'r',linewidth=3,label=r'P_1_x')
    lns_particle_1_pos_y_t = ax_particle_positions_time.plot(time, p1_data_mat[:,2],'g',linewidth=3,label=r'P_1_y')
    lns_particle_1_pos_z_t = ax_particle_positions_time.plot(time, p1_data_mat[:,3],'b',linewidth=3,label=r'P_1_z')
    ax_particle_positions_time.set_title("Position of particle")
    ax_particle_positions_time.set_xlabel('Time [s]')
    ax_particle_positions_time.set_ylabel('Position [m]')
    ax_particle_positions_time.legend()
    figure_partcle_positions_time.tight_layout()

    # ==PARTICLE POSITIONS AND RADII IN X=======================================================================================
    figure_particle_positions_and_radii, ax_particle_positions_and_radii = plt.subplots()
    lns_particle_1_pos = ax_particle_positions_and_radii.plot(time, p1_data_mat[:,1], 'r*', label=r'P_1_x')
    lns_particle_1_radii = ax_particle_positions_and_radii.plot(time, p1_data_mat[:,1]+p1_data_mat[:,7], 'r')
    lns_particle_2_pos = ax_particle_positions_and_radii.plot(time, p2_data_mat[:,1], 'b*', label=r'P_2_x')
    lns_particle_2_radii = ax_particle_positions_and_radii.plot(time, p2_data_mat[:,1]-p2_data_mat[:,7], 'b')
    ax_particle_positions_and_radii.set_title('Position and radii of particles')
    ax_particle_positions_and_radii.set_xlabel('Time [s]')
    ax_particle_positions_and_radii.set_ylabel('X Position [m]')
    ax_particle_positions_and_radii.legend()
    figure_particle_positions_and_radii.tight_layout()

# ==PARTICLE FORCE ALL DIRECTIONS=======================================================================================
    fig_all_particle_forces, ax_all_particle_forces = plt.subplots()
    lns_particle_1_force_x = ax_all_particle_forces.plot(time, p1_data_mat[:, 10], 'r-', linewidth=3, label=r'F_x')
    lns_particle_1_force_y = ax_all_particle_forces.plot(time, p1_data_mat[:, 11], 'g-', linewidth=3, label=r'F_y')
    lns_particle_1_force_z = ax_all_particle_forces.plot(time, p1_data_mat[:, 12], 'b-', linewidth=3, label=r'F_z')
    lns_particle_1_force_res = ax_all_particle_forces.plot(
        time, (p1_data_mat[:, 10]**2 + p1_data_mat[:, 11]**2 + p1_data_mat[:, 12]**2)**(1/2),
        'c--', linewidth=3, label=r'F_res')
    ax_all_particle_forces.set_title("Forces on particles")
    ax_all_particle_forces.set_xlabel("Simulation time [s]")
    ax_all_particle_forces.set_ylabel("Force on particle [N]")
    ax_all_particle_forces.legend()
    fig_all_particle_forces.tight_layout()

    # ==CONTACT FORCE ALL DIRECTIONS====================================================================================
    fig_all_contact_forces, ax_all_contact_forces = plt.subplots()
    lns_contact_1_force_n = ax_all_contact_forces.plot(time, contact_data_mat[:, 6], 'y-', linewidth=3, label=r'F_n')
    lns_contact_1_force_x = ax_all_contact_forces.plot(time, contact_data_mat[:, 10], 'r-', linewidth=3, label=r'F_Tx')
    lns_contact_1_force_y = ax_all_contact_forces.plot(time, contact_data_mat[:, 11], 'g-', linewidth=3, label=r'F_Ty')
    lns_contact_1_force_z = ax_all_contact_forces.plot(time, contact_data_mat[:, 12], 'b-', linewidth=3, label=r'F_Tz')
    lns_contact_1_force_res = ax_all_contact_forces.plot(
        time, (contact_data_mat[:, 10] ** 2 + contact_data_mat[:, 11] ** 2 + contact_data_mat[:, 12] ** 2) ** (1 / 2),
        'c-', linewidth=3, label=r'F_Tres')
    ax_all_contact_forces.set_title("Force in contact")
    ax_all_contact_forces.set_xlabel("Simulation time [s]")
    ax_all_contact_forces.set_ylabel("Force in particle [N]")
    ax_all_contact_forces.legend()
    fig_all_contact_forces.tight_layout()

    # ==FORCE DISPLACEMENT==============================================================================================
    fig_force_disp, ax_force_disp = plt.subplots()
    lns_p1_force_disp = ax_force_disp.plot(
        p1_data_mat[:,1]-p1_data_mat[0,1], -p1_data_mat[:, 10], 'y-', linewidth=3, label=r'F_x')
    # lns_contact_1_force_x = ax_all_contact_forces.plot(time, contact_data_mat[:, 10], 'r-', linewidth=3, label=r'F_Tx')
    # lns_contact_1_force_y = ax_all_contact_forces.plot(time, contact_data_mat[:, 11], 'g-', linewidth=3, label=r'F_Ty')
    # lns_contact_1_force_z = ax_all_contact_forces.plot(time, contact_data_mat[:, 12], 'b-', linewidth=3, label=r'F_Tz')
    # lns_contact_1_force_res = ax_all_contact_forces.plot(
    #     time, (contact_data_mat[:, 10] ** 2 + contact_data_mat[:, 11] ** 2 + contact_data_mat[:, 12] ** 2) ** (1 / 2),
    #     'c-', linewidth=3, label=r'F_Tres')
    ax_force_disp.set_title("Force in particle")
    ax_force_disp.set_xlabel("Displacement [m]")
    ax_force_disp.set_ylabel("Force in particle [N]")
    ax_force_disp.legend()
    # ax_force_disp.tight_layout()

    print('Showing Plot')
    plt.show()