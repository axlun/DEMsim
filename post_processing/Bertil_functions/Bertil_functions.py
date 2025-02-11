import paramiko
import sys

sys.path.insert(0, 'c:/Users/Axel/Documents/Secrets')
from secrets import pw

import numpy as np


# Input file dir and name a
# Output Data of file with each row in list
def one_file_reader(dir):
    client = paramiko.SSHClient()
    client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    client.connect("bertil.hallf.kth.se", username="axlun", password=pw)
    sftp_client = client.open_sftp()

    try:
        file = sftp_client.open(dir)
    except IOError:
        print('No such file')
        return IOError
    else:
        file_data = file.readlines()
        file.close()
        client.close()
        return file_data

    # file = sftp_client.open(dir)
    # file_data = file.readlines()
    # file.close()
    # client.close()


# Input: String with input to bertil teminal
# Output: Prints command, error string and output string in terminal
def commad_input(command0):
    print('input is:\n' + command0)
    temp = input('Enter y to send input')
    if temp.startswith('y') or temp.startswith('Y'):
        print('Sending input to Bertil')
        client = paramiko.SSHClient()
        client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        client.connect("bertil.hallf.kth.se", username="axlun", password=pw)
        stdin0, stdout0, stderr0 = client.exec_command(command0)
        print(command0)
        print(stderr0.readlines())
        print(stdout0.readlines())
        client.close()
    else:
        print('Input not sent')
    return


# Input: String with input to bertil teminal
# Output: String printed in bertil terminal
def input_output_return(command0):
    print('Sending input to Bertil:' + str(command0))
    client = paramiko.SSHClient()
    client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    client.connect("bertil.hallf.kth.se", username="axlun", password=pw)
    stdin0, stdout0, stderr0 = client.exec_command(command0)
    bertil_output = stdout0.readlines()
    client.close()
    return bertil_output

# Input: String with directory of simulation
# Output: Time, binder contacts and active binder contacts and fractured binder contacts as np.array
def fractured_binder_counter(sim_dir):
    bertil_post_process_call = 'cd /scratch/users/axlun/DEMsim/post_processing/Bertil_functions/' \
                               ' \npython3 binder_fracture_processing.py ' + sim_dir + '/contacts'
    print('Calling post processing script on Bertil')
    bertil_return = input_output_return(bertil_post_process_call)
    bertil_return = bertil_return[0].split('] [')
    time_vec = bertil_return[0][1:].split(', ')
    binder_contact_vec = bertil_return[1].split(', ')
    active_binder_contact_vec = bertil_return[2].split(', ')
    fractured_binder_contact_vec = bertil_return[3][:-2].split(', ')
    time_vec = np.array(time_vec, dtype=float)
    binder_contact_vec = np.array(binder_contact_vec, dtype=float)
    active_binder_contact_vec = np.array(active_binder_contact_vec, dtype=float)
    fractured_binder_contact_vec = np.array(fractured_binder_contact_vec, dtype=float)
    return time_vec, binder_contact_vec, active_binder_contact_vec, fractured_binder_contact_vec


# Input: String with directory of simulation
# Output: Time, particle contacts and binder contacts as np.array
def contact_counter_bertil(sim_dir):
    bertil_post_process_call = 'cd /scratch/users/axlun/DEMsim/post_processing/Bertil_functions/' \
                               ' \npython3 contact_processing.py ' + sim_dir + '/contacts'
    print('Calling post processing script on Bertil')
    bertil_return = input_output_return(bertil_post_process_call)
    bertil_return = bertil_return[0].split('] [')
    time_vec = bertil_return[0][1:].split(', ')
    binder_contact_vec = bertil_return[1].split(', ')
    particle_contact_vec = bertil_return[2].split(', ')
    binder_particle_contact_vec = bertil_return[3][:-2].split(', ')
    time_vec = np.array(time_vec, dtype=float)
    particle_contact_vec = np.array(particle_contact_vec, dtype=float)
    binder_contact_vec = np.array(binder_contact_vec, dtype=float)
    binder_particle_contact_vec = np.array(binder_particle_contact_vec, dtype=float)
    return time_vec, particle_contact_vec, binder_contact_vec, binder_particle_contact_vec


# Input: String with directory of simulation
# Output: Time as np.array and particle contact valency and binder contact valency as np.ndarray
def contact_valency_bertil(sim_dir):
    bertil_post_process_call = 'cd /scratch/users/axlun/DEMsim/post_processing/Bertil_functions/' \
                               ' \npython3 contact_valency_processing.py ' + sim_dir + '/contacts'
    print('Calling Bertil post processing script')
    bertil_return = input_output_return(bertil_post_process_call)
    time_vec = bertil_return[0][1:-2].split(', ')
    time_vec = np.array(time_vec, dtype=float)
    contact_valency_vec = bertil_return[1:]
    particle_valency_vec = contact_valency_vec[:(len(contact_valency_vec) // 2)]
    binder_valency_vec = contact_valency_vec[(len(contact_valency_vec) // 2):]

    for x in range(0, len(contact_valency_vec) // 2):
        binder_valency_vec[x] = binder_valency_vec[x][2:-2]
        particle_valency_vec[x] = particle_valency_vec[x][2:-2]
    particle_valency_vec[-1] = particle_valency_vec[-1][:-1]
    binder_valency_vec[-1] = binder_valency_vec[-1][:-1]

    particle_contact_valency_mat = []
    binder_contact_valency_mat = []
    for x in range(0, len(contact_valency_vec) // 2):
        particle_contact_valency_mat.append(particle_valency_vec[x].split())
        binder_contact_valency_mat.append(binder_valency_vec[x].split())
    return time_vec, particle_contact_valency_mat, binder_contact_valency_mat


def sum_normal_tangential_force_bertil(sim_dir):
    bertil_post_process_call = 'cd /scratch/users/axlun/DEMsim/post_processing/Bertil_functions/' \
                               ' \npython3 sum_normal_tangential_force_processing.py ' + sim_dir + '/contacts'
    print('Calling Bertil post processing script')
    bertil_return = input_output_return(bertil_post_process_call)
    bertil_return = bertil_return[0].split('] [')
    time_vec = bertil_return[0][1:].split(', ')
    time_vec = np.array(time_vec, dtype=float)

    normal_force_vec = bertil_return[1].split(', ')
    normal_force_vec = np.array(normal_force_vec, dtype=float)

    compressive_normal_force_vec = bertil_return[2].split(', ')
    compressive_normal_force_vec = np.array(compressive_normal_force_vec, dtype=float)

    tensile_normal_force_vec = bertil_return[3].split(', ')
    tensile_normal_force_vec = np.array(tensile_normal_force_vec, dtype=float)

    binder_force_vec = bertil_return[4].split(', ')
    binder_force_vec = np.array(binder_force_vec, dtype=float)

    compressive_binder_force_vec = bertil_return[5].split(', ')
    compressive_binder_force_vec = np.array(compressive_binder_force_vec, dtype=float)

    tensile_binder_force_vec = bertil_return[6].split(', ')
    tensile_binder_force_vec = np.array(tensile_binder_force_vec, dtype=float)

    particle_force_vec = bertil_return[7].split(', ')
    particle_force_vec = np.array(particle_force_vec, dtype=float)

    tangential_force_vec = bertil_return[8][:-2].split(', ')
    tangential_force_vec = np.array(tangential_force_vec, dtype=float)
    return time_vec, normal_force_vec, compressive_normal_force_vec, tensile_normal_force_vec, binder_force_vec, \
        compressive_binder_force_vec, tensile_binder_force_vec, particle_force_vec, tangential_force_vec


def bertil_contact_tracker(sim_dir, p1, p2):
    bertil_post_process_call = 'cd /scratch/users/axlun/DEMsim/post_processing/Bertil_functions/ ' \
                               '\npython3 contact_tracker_processing.py '\
                               + sim_dir + '/contacts ' + str(p1) + ' ' + str(p2)
    print('Calling Bertil post processing script')
    bertil_return = input_output_return(bertil_post_process_call)
    time_vec = bertil_return[0][1:-2].split(', ')
    bertil_return = bertil_return[1:]
    time_vec = np.array(time_vec, dtype=float)
    contact_data_vec = []

    for i in range(0, len(time_vec)):
        step_vec = bertil_return[21 * i][:-1]
        for j in range(1, 21):
            step_vec += ' ' + bertil_return[21 * i + j][:-1]
        # print(step_vec)
        if i == 0:
            contact_data_vec = np.fromstring(step_vec, dtype=float, sep=' ')
            continue
        contact_data_vec = np.vstack((contact_data_vec, np.fromstring(step_vec, dtype=float, sep=' ')))

    return time_vec, contact_data_vec

def bertil_single_contact_tracker(sim_dir, p1):
    bertil_post_process_call = 'cd /scratch/users/axlun/DEMsim/post_processing/Bertil_functions/ ' \
                               '\npython3 single_contact_tracker_processing.py '\
                               + sim_dir + '/contacts ' + str(p1)
    print('Calling Bertil post processing script')
    bertil_return = input_output_return(bertil_post_process_call)
    # print(bertil_return)
    time_vec = bertil_return[0][1:-2].split(', ')
    time_vec = np.array(time_vec, dtype=float)
    contact_vec = bertil_return[1][1:-2].split(',')
    contact_vec = np.array(contact_vec,dtype=float)
    columns = 2
    contact_vec= contact_vec.reshape((len(contact_vec)//columns,columns))

    # bertil_return = bertil_return[1:]
    # contact_data_vec = []

    # for i in range(0, len(time_vec)):
    #     step_vec = bertil_return[3 * i][:-1]
    #     for j in range(1, 3):
    #         step_vec += ' ' + bertil_return[3 * i + j][:-1]
    #     # print(step_vec)
    #     if i == 0:
    #         contact_data_vec = np.fromstring(step_vec, dtype=float, sep=' ')
    #         continue
    #     contact_data_vec = np.vstack((contact_data_vec, np.fromstring(step_vec, dtype=float, sep=' ')))

    return time_vec, contact_vec

def bertil_single_particle_tracker(sim_dir, p1):
    bertil_post_process_call = 'cd /scratch/users/axlun/DEMsim/post_processing/Bertil_functions/ ' \
                               '\npython3 single_particle_tracker_processing.py '\
                               + sim_dir + '/particles ' + str(p1)
    print('Calling Bertil post processing script')
    bertil_return = input_output_return(bertil_post_process_call)
    print(bertil_return)
    print(len(bertil_return))
    time_vec = bertil_return[0][1:-2].split(', ')
    time_vec = np.array(time_vec, dtype=float)

    particle_vec = bertil_return[1][1:-2].split(',')
    particle_vec = np.array(particle_vec,dtype=float)
    columns = 13
    particle_vec= particle_vec.reshape((len(particle_vec)//columns,columns))
    return time_vec, particle_vec


def bertil_layer_height(sim_dir, no_particle_in_avg=10):
    bertil_post_process_call = 'cd /scratch/users/axlun/DEMsim/post_processing/Bertil_functions/ ' \
                               '\npython3 layer_height_processing.py ' \
                                + sim_dir + '/particles ' + str(no_particle_in_avg)
    print('Calling Bertil layer height post processing script')
    bertil_return = input_output_return(bertil_post_process_call)
    time_vec = bertil_return[0][1:-2].split(', ')
    time_vec = np.array(time_vec, dtype=float)

    particle_vec = bertil_return[1][1:-2].split(',')
    particle_vec = np.array(particle_vec,dtype=float)
    return time_vec, particle_vec

# Input: String with directory of simulation
# Output: np.array of force_data, surface_force_index, surface_position_index,surface_position_data, periodic_BC_data, force_fabric_tensor_data as floats
def bertil_data_gatherer(simulation_directory):
    #================SURFACE FORCE DATA================================================================================
    print('Getting surface force data from bertil')

    force_data_list = one_file_reader(simulation_directory + '/surface_forces.dou')
    if force_data_list == IOError:
        force_data_list = []
    else:
        force_data_list = np.asarray(force_data_list)
        force_data_list = np.asarray(force_data_list)
        force_data_list = np.char.split(force_data_list, ', ')
        first_line = force_data_list[0]
        id_idx = [i for i in range(len(first_line)) if first_line[i].startswith('ID')]
        surface_force_index = id_idx
        for i in force_data_list[:]:
            if 'force_data' in locals():
                force_data = np.vstack((force_data, np.asarray(i)))
            else:
                force_data = np.asarray(i)
        for i in enumerate(force_data): force_data[i[0], -1] = force_data[i[0], -1][:-1]
    #================PERIODIC BC DATA===================================================================================
    print('Getting periodic BC data from bertil')

    periodic_BC_data_list = one_file_reader(simulation_directory + '/periodic_bc.dou')
    if periodic_BC_data_list == IOError:
        periodic_BC_data = []
    else:
        periodic_BC_data_list = np.asarray(periodic_BC_data_list)
        periodic_BC_data_list = np.char.split(periodic_BC_data_list, ', ')
        for i in periodic_BC_data_list[:]:
            if 'periodic_BC_data' in locals():
                periodic_BC_data = np.vstack((periodic_BC_data, np.asarray(i)))
            else:
                periodic_BC_data = np.asarray(i)
        for i in enumerate(periodic_BC_data): periodic_BC_data[i[0], -1] = periodic_BC_data[i[0], -1][:-1]
        periodic_BC_data = periodic_BC_data.astype(float)

    # ================SURFACE POSITION DATA=============================================================================
    print('Getting surface position data from bertil')
    surface_position_data_list = one_file_reader(simulation_directory + '/surface_positions.dou')

    first_line = surface_position_data_list[0].split(', ')
    id_idx = [i for i in range(len(first_line)) if first_line[i].startswith('ID')]
    surface_position_index = id_idx
    surface_types = [first_line[idx + 1][5:] for idx in id_idx]
    surface_indices = [i for i, surface_type in enumerate(surface_types) if surface_type == 'PointSurface']

    surface_position_data_list = np.asarray(surface_position_data_list)
    surface_position_data_list = np.char.split(surface_position_data_list, ', ')
    for i in surface_position_data_list[:]:
        if 'surface_position_data' in locals():
            surface_position_data = np.vstack((surface_position_data, np.asarray(i)))
        else:
            surface_position_data = np.asarray(i)
    for i in enumerate(surface_position_data):
        surface_position_data[i[0], -1] = surface_position_data[i[0], -1][:-1]
    # ================FORCE FABRIC TENSOR DATA==========================================================================
    print('Getting force fabric tensor data from bertil')
    force_fabric_tensor_data_list = one_file_reader(simulation_directory + '/force_fabric_tensor.dou')
    force_fabric_tensor_data_list = np.asarray(force_fabric_tensor_data_list)
    force_fabric_tensor_data_list = np.char.split(force_fabric_tensor_data_list, ', ')
    for i in force_fabric_tensor_data_list[:]:
        if 'force_fabric_tensor_data' in locals():
            force_fabric_tensor_data = np.vstack((force_fabric_tensor_data, np.asarray(i)))
        else:
            force_fabric_tensor_data = np.asarray(i)
    for i in enumerate(force_fabric_tensor_data): force_fabric_tensor_data[i[0], -1] = force_fabric_tensor_data[
                                                                                           i[0], -1][:-1]
    force_fabric_tensor_data = force_fabric_tensor_data.astype(float)

    # ====================KINETIC ENERGY================================================================================
    print('Getting kinetic energy data from bertil')
    kinetic_energy_data_list = one_file_reader(simulation_directory + '/kinetic_energy.dou')
    kinetic_energy_data_list = np.asarray(kinetic_energy_data_list)
    kinetic_energy_data_list = np.char.split(kinetic_energy_data_list, ', ')
    for i in kinetic_energy_data_list[:]:
        if 'kinetic_energy_data' in locals():
            kinetic_energy_data = np.vstack((kinetic_energy_data, np.asarray(i)))
        else:
            kinetic_energy_data = np.asarray(i)
    for i in enumerate(kinetic_energy_data): kinetic_energy_data[i[0], -1] = kinetic_energy_data[i[0], -1][:-1]
    kinetic_energy_data = kinetic_energy_data.astype(float)

    return force_data, surface_force_index, surface_position_index, surface_position_data, periodic_BC_data, force_fabric_tensor_data, kinetic_energy_data

def fractured_particle_gatherer(simulation_directory):

    print("Getting fractured particles from Bertil.")
    fractured_particle_list = one_file_reader(simulation_directory + '/fractured_particles.dou')
    # fractured_particle_list = one_file_reader(simulation_directory + '/kinetic_energy.dou')

    data = []
    time = np.array([])
    for line in fractured_particle_list:
        row = line.strip().split(',')
        row = np.array([float(x) for x in row])  # Create a NumPy array for each row
        data.append(row)
        time= np.append(time, row[-1])
    fracture_array = np.array(data, dtype=object)

    number_of_fractures = np.array([np.shape(x)[0]-1 for x in fracture_array])
    return time, number_of_fractures, fracture_array


if __name__ == '__main__':
    command0 = "cd /scratch/users/axlun/DEMsim/"
    command = "ls"

    input_command = command0 + " \n " + command
    commad_input(input_command)
