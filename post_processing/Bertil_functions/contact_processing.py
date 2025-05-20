import sys
import os
import re
import pandas as pd

if __name__ == '__main__':
    argument_string = sys.argv[1]

    contact_files = os.listdir(argument_string)
    time = []
    contact_time_and_file_name_dict = {}
    for i in range(0, len(contact_files)):
        contact_file = contact_files[i]
        time_stamp = re.split(r'\Acontacts_', contact_file)
        time_stamp = re.split(r'.dou\Z', time_stamp[1])

        time.append(float(time_stamp[0]))
        contact_time_and_file_name_dict[time_stamp[0]] = contact_file

    time.sort()
    time_vec = []
    binder_contact_vec = []
    particle_contact_vec = []
    binder_particle_contact_vec = []

    for i in range(0, len(time)):
        binder_contact = 0
        particle_contact = 0
        key = str(time[i])
        if time[i].is_integer():
            key = str(int(time[i]))
        ## Read file here
        file_to_open = argument_string + '/' + contact_time_and_file_name_dict[key]
        data = pd.read_csv(file_to_open, header=None).to_numpy()
        binder_contact = 0
        particle_contact = 0
        binder_particle_contact = 0
        for j in range(0, len(data[:, 0])):
            if data[j, 6] != 0 and (data[j, 6]) == (data[j, 7]):
                binder_contact += 1
            if data[j, 6] != 0 and (data[j, 6]) == (data[j, 8]):
                particle_contact += 1
            # if data[j, 6] != 0 and data[j, 7] != 0 and data[j, 8] != 0 and data[j, 6] == data[j, 7] + data[j, 8]:
            if data[j, 6] != 0 and data[j, 7] != 0 and data[j, 8] != 0 and \
                    abs((data[j, 6] - (data[j, 7] + data[j, 8])) / data[j, 6]) < 1e-3:
                binder_particle_contact += 1
                # ==USE IF COMBINED CONTACT SHOULD CONTRIBUTE TO TOTAL NUMBER OF CONTACTS===================
                binder_contact += 1
                particle_contact += 1

        binder_contact_vec.append(binder_contact)
        particle_contact_vec.append(particle_contact)
        binder_particle_contact_vec.append(binder_particle_contact)
        time_vec.append(float(key))

    print(time_vec, binder_contact_vec, particle_contact_vec, binder_particle_contact_vec)
