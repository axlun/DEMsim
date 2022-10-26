import sys
import os
import re
import numpy as np
import pandas as pd
np.set_printoptions(threshold=np.inf)

if __name__ == '__main__':
    argument_string = sys.argv[1]

#==================GET NUMBER OF PARTICLES========================
    particle_files = os.listdir(argument_string.replace("contacts","particles"))

    with open(argument_string.replace("contacts","particles")+'/'+particle_files[-1]) as opened_particle_file:
        lines = opened_particle_file.readlines()
    number_of_particles = int(lines[-1].split(', ')[0]) - int(lines[0].split(', ')[0])+1
    max_wall_index = int(lines[0].split(', ')[0])
    contact_files = os.listdir(argument_string)
    time = []
    contact_time_and_file_name_dict = {}
#   Get all time steps for contacts results
    for i in range(0,len(contact_files)):
        contact_file = contact_files[i]
        time_stamp = re.split(r'\Acontacts_',contact_file)
        time_stamp = re.split(r'.dou\Z', time_stamp[1])
        time.append(float(time_stamp[0]))
        contact_time_and_file_name_dict[time_stamp[0]] = contact_file

    number_of_contacts_observed = 8
    time.sort()
    binder_contact_vec = []
    particle_contact_vec = []
    particle_results_vec = np.zeros(shape=(len(time),number_of_contacts_observed),dtype=int,order='C')
    binder_results_vec = np.zeros(shape=(len(time),number_of_contacts_observed),dtype=int,order='C')
    for i in range(0,len(time)):
        particle_array = [0]* number_of_particles
        particle_binder_array = [0]* number_of_particles
        key = str(time[i])
        if time[i].is_integer():
            key = str(int(time[i]))

        ## Read file here
        file_to_open = argument_string+'/'+contact_time_and_file_name_dict[key]
        particle_contact_number_vec = [0] * number_of_contacts_observed
        binder_contact_number_vec = [0] * number_of_contacts_observed

        line_data = pd.read_csv(file_to_open).to_numpy()
        for j in range(0, len(line_data[:,0])):
            # PARTICLE CONTACT
            if float(line_data[j,5]) >= 0:
                particle_array[int(line_data[j,0]) - max_wall_index] += 1
                if int(line_data[j,1]) >= max_wall_index:
                    particle_array[int(line_data[j,1]) - max_wall_index] += 1
            #BINDER CONTACT
            if float(line_data[j,5]) < 0 and float(line_data[j,6]) != 0.:
                particle_binder_array[int(line_data[j,0]) - max_wall_index] += 1
                if int(line_data[j,1]) >= max_wall_index:
                    particle_binder_array[int(line_data[j,1]) - max_wall_index] += 1

        #=OLD READING METHOD===========================================================================================
        # with open(file_to_open) as opened_contact_file:
        #     lines = opened_contact_file.readlines()
        # for j in range(0, len(lines)):
        #     line_data = lines[j].split(', ')
        #     if float(line_data[5]) >= 0:
        #         particle_array[int(line_data[0]) - max_wall_index] += 1
        #         if int(line_data[1]) >= max_wall_index:
        #             particle_array[int(line_data[1]) - max_wall_index] += 1
        #
        #     if float(line_data[5]) < 0 and float(line_data[6]) != 0.:
        #         particle_binder_array[int(line_data[0]) - max_wall_index] += 1
        #         if int(line_data[1]) >= max_wall_index:
        #             particle_binder_array[int(line_data[1]) - max_wall_index] += 1
        #==============================================================================================================

        for k in particle_array:
            if k > number_of_contacts_observed - 1:
                k = number_of_contacts_observed - 1
            particle_contact_number_vec[k] += 1
        for l in particle_binder_array:
            if l > number_of_contacts_observed - 1:
                l = number_of_contacts_observed - 1
            binder_contact_number_vec[l] += 1
        particle_results_vec[i, :] = particle_contact_number_vec
        binder_results_vec[i, :] = binder_contact_number_vec
    print(time)
    print(particle_results_vec)
    print(binder_results_vec)
#     print('Script end')