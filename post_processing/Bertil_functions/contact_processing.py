import sys
import os
import re


if __name__ == '__main__':
    argument_string = sys.argv[1]

    # print('Hello world')

    # print(argument_string)
    contact_files = os.listdir(argument_string)
#    print(contact_files)
    time = []
    contact_time_and_file_name_dict  = {}
    for i in range(0,len(contact_files)):
        contact_file = contact_files[i]
        time_stamp = re.split(r'\Acontacts_',contact_file)
#        print(time_stamp[1])
        time_stamp = re.split(r'.dou\Z', time_stamp[1])
#        print(time_stamp)
#        contact_file = re.split(r'\n\Z',contact_file)
#        print(contact_file)
        time.append(float(time_stamp[0]))
        contact_time_and_file_name_dict[time_stamp[0]] = contact_file

    # print(contact_time_and_file_name_dict)
    time.sort()
    # print(time)
    time_vec = []
    binder_contact_vec = []
    particle_contact_vec = []
    for i in range(0,len(time)):
        binder_contact = 0
        particle_contact = 0
        key = str(time[i])
        if time[i].is_integer():
            key = str(int(time[i]))
        # print(key)

        ## Read file here
        file_to_open = argument_string+'/'+contact_time_and_file_name_dict[key]
        with open(file_to_open) as opened_contact_file:
            lines = opened_contact_file.readlines()
        # print(lines)

        for j in range(0, len(lines)):
             line_data = lines[j].split(', ')
             if float(line_data[7]) != 0 and float(line_data[8]) == 0.0:
                # print('Binder flag')
                binder_contact += 1
             if float(line_data[8]) != 0:
                # print('Particle flag')
                particle_contact += 1

            # if float(line_data[6]) > 0:
            #     print(line_data[6])
        binder_contact_vec.append(binder_contact)
        particle_contact_vec.append(particle_contact)
        time_vec.append(float(key))
        # print('Time: ' + key)
        # print('Binder contacts: ' + str(binder_contact))
        # print('Particle contacts: ' + str(particle_contact))


    print(time_vec,binder_contact_vec,particle_contact_vec)
    # print('Script end')