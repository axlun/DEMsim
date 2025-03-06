import sys
import os
import re
import pandas as pd

if __name__ == '__main__':
    argument_string = sys.argv[1]
    contact_files = os.listdir(argument_string)
    time = []
    contact_time_and_file_name_dict = {}
    for i in range(0,len(contact_files)):
        contact_file = contact_files[i]
        time_stamp = re.split(r'\Acontacts_',contact_file)
        time_stamp = re.split(r'.dou\Z', time_stamp[1])
        time.append(float(time_stamp[0]))
        contact_time_and_file_name_dict[time_stamp[0]] = contact_file
    time.sort()
    time_vec = []
    binder_contact_vec = []
    active_binder_contact_vec = []
    fractured_binder_contact_vec = []
    for i in range(0,len(time)):
        key = str(time[i])
        if time[i].is_integer():
            key = str(int(time[i]))
        active_binder_contact = 0
        binder_contact = 0
        fractured_binder_contact = 0
        file_to_open = argument_string + '/' + contact_time_and_file_name_dict[key]
        try:
            data = pd.read_csv(file_to_open,header=None).to_numpy()
            for j in range(0, len(data[:, 0])):
                if data[j, 6] != 0 and (data[j, 6]) == (data[j, 7]):
                    active_binder_contact += 1
                if data[j, 21] == 1 :
                    fractured_binder_contact += 1
                if data[j, 19] == 1 :
                    binder_contact += 1
        except:
           data = data
        binder_contact_vec.append(binder_contact)
        active_binder_contact_vec.append(active_binder_contact)
        fractured_binder_contact_vec.append(fractured_binder_contact)
        time_vec.append(float(key))

    print(time_vec,binder_contact_vec,active_binder_contact_vec,fractured_binder_contact_vec)
