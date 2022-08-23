
import paramiko
import sys
sys.path.insert(0, 'c:/Users/Axel/Documents/Secrets')
from secrets import pw

import numpy as np

#Input file dir and name a
#Output Data of file with each row in list
def one_file_reader(dir):
    client = paramiko.SSHClient()
    client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    client.connect("bertil.hallf.kth.se", username="axlun", password=pw)
    sftp_client = client.open_sftp()
    file = sftp_client.open(dir)
    file_data = file.readlines()
    file.close()
    client.close()
    return file_data

#Input: String with input to bertil teminal
#Output: Prints command, error string and output string in terminal
def commad_input(command0):
    print('input is:\n'+command0)
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
        # stdin, stdout, stderr = client.exec_command(command)
        # print(stderr.readlines())
        # lines = stdout.readlines()
        # print(lines)
        client.close()
    else:
        print('Input not sent')
    return


#Input: String with input to bertil teminal
#Output: String printed in bertil terminal
def input_output_return(command0):
#    print('input is:\n'+command0)
    temp = 'Y'#input('Enter y to send input')
    if temp.startswith('y') or temp.startswith('Y'):
#        print('Sending input to Bertil')
        client = paramiko.SSHClient()
        client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        client.connect("bertil.hallf.kth.se", username="axlun", password=pw)

        stdin0, stdout0, stderr0 = client.exec_command(command0)
#        print(command0)
#        print(stderr0.readlines())
        bertil_output = stdout0.readlines()
#        print(bertil_output)
        # stdin, stdout, stderr = client.exec_command(command)
        # print(stderr0.readlines())
        # lines = stdout.readlines()
        # print(lines)
        client.close()
    else:
        print('Input not sent')
        return
#    print(bertil_output)
    return bertil_output


#Input: String with directory of simulation
#Output: Time, particle contacts and binder contacts as np.array
def contact_counter_bertil(sim_dir):
    bertil_post_process_call = 'cd /scratch/users/axlun/DEMsim/post_processing/Bertil_functions/ \npython3 contact_processing.py ' + sim_dir+'/contacts'
    bertil_return = input_output_return(bertil_post_process_call)
    bertil_return = bertil_return[0].split('] [')
    time_vec = bertil_return[0][1:].split(', ')
    binder_contact_vec = bertil_return[1].split(', ')
    particle_contact_vec = bertil_return[2][:-2].split(', ')
    time_vec = np.array(time_vec, dtype=float)
    particle_contact_vec = np.array(particle_contact_vec, dtype=float)
    binder_contact_vec = np.array(binder_contact_vec, dtype=float)
    return time_vec,particle_contact_vec,binder_contact_vec



#Input: String with directory of simulation
#Output: Time as np.array and particle contact valency and binder contact valency as np.ndarray
def contact_valency_bertil(sim_dir):
    bertil_post_process_call = 'cd /scratch/users/axlun/DEMsim/post_processing/Bertil_functions/ \npython3 contact_valency_processing.py ' + sim_dir+'/contacts'
    bertil_return = input_output_return(bertil_post_process_call)
    time_vec = bertil_return[0][1:-2].split(', ')
    time_vec = np.array(time_vec, dtype=float)
    contact_valency_vec = bertil_return[1:]
    particle_valency_vec = contact_valency_vec[:(len(contact_valency_vec)//2)]
    binder_valency_vec = contact_valency_vec[(len(contact_valency_vec)//2):]

    for x in range(0, len(contact_valency_vec)//2):
        binder_valency_vec[x] = binder_valency_vec[x][2:-2]
        particle_valency_vec[x] = particle_valency_vec[x][2:-2]
    particle_valency_vec[-1] = particle_valency_vec[-1][:-1]
    binder_valency_vec[-1] = binder_valency_vec[-1][:-1]

    particle_contact_valency_mat = []
    binder_contact_valency_mat = []
    for x in range(0,len(contact_valency_vec)//2):
        particle_contact_valency_mat.append(particle_valency_vec[x].split())
        binder_contact_valency_mat.append(binder_valency_vec[x].split())
    return time_vec,particle_contact_valency_mat,binder_contact_valency_mat


#Input: String with directory of simulation
#Output: np.array of force_data, surface_force_index, surface_position_index,surface_position_data, periodic_BC_data, force_fabric_tensor_data as floats
def bertil_data_gatherer(simulation_directory):
    ##================SURFACE FORCE DATA===================================================================================
    force_data_list = one_file_reader(simulation_directory + '/surface_forces.dou')
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
    ##================PERIODIC BC DATA==================================================================================
    periodic_BC_data_list = one_file_reader(simulation_directory + '/periodic_bc.dou')
    periodic_BC_data_list = np.asarray(periodic_BC_data_list)
    periodic_BC_data_list = np.char.split(periodic_BC_data_list, ', ')
    for i in periodic_BC_data_list[:]:
        if 'periodic_BC_data' in locals():
            periodic_BC_data = np.vstack((periodic_BC_data, np.asarray(i)))
        else:
            periodic_BC_data = np.asarray(i)
    for i in enumerate(periodic_BC_data): periodic_BC_data[i[0], -1] = periodic_BC_data[i[0], -1][:-1]
    periodic_BC_data = periodic_BC_data.astype(float)

    # ================SURFACE POSITION DATA==========================================================================
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
    ##================FORCE FABRIC TENSOR DATA==================================================================================
    force_fabric_tensor_data_list = one_file_reader(simulation_directory + '/force_fabric_tensor.dou')
    force_fabric_tensor_data_list = np.asarray(force_fabric_tensor_data_list)
    force_fabric_tensor_data_list = np.char.split(force_fabric_tensor_data_list, ', ')
    for i in force_fabric_tensor_data_list[:]:
        if 'force_fabric_tensor_data' in locals():
            force_fabric_tensor_data = np.vstack((force_fabric_tensor_data, np.asarray(i)))
        else:
            force_fabric_tensor_data = np.asarray(i)
    for i in enumerate(force_fabric_tensor_data): force_fabric_tensor_data[i[0], -1] = force_fabric_tensor_data[i[0], -1][:-1]
    force_fabric_tensor_data = force_fabric_tensor_data.astype(float)



    return force_data, surface_force_index, surface_position_index,surface_position_data, periodic_BC_data, force_fabric_tensor_data






if __name__ == '__main__':
    command0 = "cd /scratch/users/axlun/DEMsim/"
    command = "ls"

    input_command = command0 +" \n " + command
    commad_input(input_command)
