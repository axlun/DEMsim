import os

import numpy as np
import re
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.style
import pandas as pd


def get_parameter(parameter_dir, parameter):
    parameter_array = np.genfromtxt(parameter_dir, dtype='U30,f', delimiter='=')
    for index, item in np.ndenumerate(parameter_array):
        if parameter_array[index]['f0'] == parameter:
            return parameter_array[index]['f1']


if __name__ == '__main__':
    plt.style.use('axel_style')

    # simulation_directory = "C:/Users/Axel/Documents/DEM/results/particle_tests/fracture_particle_plate_compression/SN_1/"
    # simulation_directory = "C:/Users/Axel/Documents/DEM/results/particle_tests/fracture_particle_plate_compression/SN_2/"
    # simulation_directory = "C:/Users/Axel/Documents/DEM/results/particle_tests/fracture_particle_plate_compression/SN_3/"
    simulation_directory = "C:/Users/Axel/Documents/DEM/results/particle_tests/fracture_particle_plate_compression/SN_4/"
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
    contact_1_data = []
    contact_2_data = []
    for i in range(0,len(time)):
        key = str(time[i])
        if time[i].is_integer():
            key = str(int(time[i]))
        data = np.genfromtxt(simulation_directory+'particles/'+particle_time_and_file_name_dict[key],delimiter=', ')
        p1_data.append(data)
        contact_data = np.genfromtxt((simulation_directory+'contacts/'+particle_time_and_file_name_dict[key]).replace(
                "particles","contacts"),delimiter=', ')
        contact_1_data.append(contact_data[0, :])
        contact_2_data.append(contact_data[1, :])
    p1_data_mat = np.stack(p1_data,axis=0)
    contact_1_data_mat = np.stack(contact_1_data, axis=0)
    contact_2_data_mat = np.stack(contact_2_data, axis=0)
    force_fabric_tensor = np.genfromtxt(simulation_directory+'force_fabric_tensor.dou',delimiter=', ')
    kinetic_energy = np.genfromtxt(simulation_directory+'kinetic_energy.dou',delimiter=', ')
    kinetic_energy_particle_1 = p1_data_mat[:,7]
    potential_energy_particle_1 = 4/3 * 3.14* p1_data_mat[:,7] ** 3 * 4800 * 10 * (p1_data_mat[:,3] - p1_data_mat[:,7])
    particle_diameter = p1_data_mat[:,7] * 2

    surface_forces_data = np.genfromtxt(simulation_directory+'surface_forces.dou', delimiter=', ')
    surface_1_force = surface_forces_data[:,1]
    surface_2_force = surface_forces_data[:,6]
    surface_positions_data = np.genfromtxt(simulation_directory+'surface_positions.dou', delimiter=', ')
    surface_1_position = surface_positions_data[:,5]

    fig_particle_forces, ax_particle_forces = plt.subplots()
    lns_particle_1_force = ax_particle_forces.plot(time,p1_data_mat[:,12],'r*-',linewidth=3,label=r'P_1')
    ax_particle_forces.set_xlabel("Simulation time [s]")
    ax_particle_forces.set_ylabel("Z Force in particle [N]")
    ax_particle_forces.legend()

    figure_partcle_positions_time, ax_particle_positions_time = plt.subplots()
    lns_particle_1_pos_t = ax_particle_positions_time.plot(time, p1_data_mat[:,3],'r',linewidth=3,label=r'P_1')
    ax_particle_positions_time.set_xlabel('Time [s]')
    ax_particle_positions_time.set_ylabel('z position [m]')
    ax_particle_positions_time.legend()

    # =====================================OVERLAP TO CONTACT FORCE=====================================================
    figure_overlap_force, ax_overlap_force = plt.subplots()
    lns_overlap_force = ax_overlap_force.plot(
        contact_1_data_mat[:,5],contact_1_data_mat[:,6],'r',linewidth=3,label=r'Contact force')
    lns_binder_force = ax_overlap_force.plot(
        contact_1_data_mat[:,5],contact_1_data_mat[:,7],'g--',linewidth=3,label=r'Binder force')
    lns_particle_force = ax_overlap_force.plot(
        contact_1_data_mat[:,5],contact_1_data_mat[:,8],'b--*',linewidth=3,label=r'Particle force')
    ax_overlap_force.set_xlabel("Overlap [m]")
    ax_overlap_force.set_ylabel("Contact force [N]")
    ax_overlap_force.legend()

    # =====================================OVERLAP TO CONTACT FORCE=====================================================
    figure_overlap_force_2, ax_overlap_force_2 = plt.subplots()
    lns_overlap_force_2 = ax_overlap_force_2.plot(
        contact_2_data_mat[:,5],contact_2_data_mat[:,6],'r',linewidth=3,label=r'Contact force')
    lns_binder_force_2 = ax_overlap_force_2.plot(
       contact_2_data_mat[:,5],contact_2_data_mat[:,7],'g--',linewidth=3,label=r'Binder force')
    lns_particle_force_2 = ax_overlap_force_2.plot(
        contact_2_data_mat[:,5],contact_2_data_mat[:,8],'b--*',linewidth=3,label=r'Particle force')
    ax_overlap_force_2.set_xlabel("Overlap [m]")
    ax_overlap_force_2.set_ylabel("Contact force [N]")
    ax_overlap_force_2.legend()

    # =====================================OVERLAP TIME=================================================================
    figure_overlap_time, ax_overlap_time = plt.subplots()
    lns_overlap_time = ax_overlap_time.plot(time, contact_1_data_mat[:, 5], 'r', linewidth=3,
                                              label=r'overlap')
    # lns_binder_force = ax_overlap_force.plot(contact_1_data_mat[:, 5], contact_1_data_mat[:, 7], 'g--', linewidth=3,
    #                                         label=r'Binder force')
    # lns_particle_force = ax_overlap_force.plot(contact_1_data_mat[:, 5], contact_1_data_mat[:, 8], 'b--', linewidth=3,
    #                                           label=r'Particle force')
    ax_overlap_time.set_ylabel("Overlap [m]")
    ax_overlap_time.set_xlabel("Time [s]")
    ax_overlap_time.legend()

    figure_disp_force, ax_disp_force = plt.subplots()
    lns_disp_force = ax_disp_force.plot((surface_1_position[0]-surface_1_position)*1E6, surface_1_force * 1E3)
    ax_disp_force.set_ylabel('Force [mN]')
    ax_disp_force.set_xlabel('Displacement [µm]')
    figure_disp_force.tight_layout()

    # EXPERIMENTAL RESULTS
    compression_experiments_df = pd.read_csv(
        "G:/My Drive/Skola/KTH/PhD/Litteratur/Fracture Testing of Lithium-Ion Battery Cathode Secondary Particles "
        "in-situ inside the Scanning Electron Microscope/flat_platen_particle_compression_load_displacement.csv")
    experiment_displacement = compression_experiments_df['displacement [µm]'].to_numpy()
    experiment_load = compression_experiments_df['load [mN]'].to_numpy()

    figure_exp_sim, ax_exp_sim = plt.subplots()
    lns_exp = ax_exp_sim.plot(experiment_displacement, experiment_load, label='exp')
    lns_sim = ax_exp_sim.plot((surface_1_position[0]-surface_1_position)*1E6, surface_1_force * 1E3, label='sim')
    ax_exp_sim.set_ylabel('Force [mN]')
    ax_exp_sim.set_xlabel('Displacement [µm]')
    ax_exp_sim.legend()
    figure_exp_sim.tight_layout()

    figure_normalised_exp_sim, ax_normalised_exp_sim = plt.subplots()
    lns_exp = ax_normalised_exp_sim.plot(experiment_displacement / 16 , experiment_load * 1E-3/((16E-6)**2) , label='exp')
    lns_sim = ax_normalised_exp_sim.plot((surface_1_position[0]-surface_1_position)/particle_diameter, surface_1_force/(particle_diameter**2), label='sim')
    ax_normalised_exp_sim.set_ylabel('F / $R^2$  [N/$m^2$]')
    ax_normalised_exp_sim.set_xlabel('u/ R [m/m]')
    ax_normalised_exp_sim.legend()
    figure_normalised_exp_sim.tight_layout()

    plt.show()