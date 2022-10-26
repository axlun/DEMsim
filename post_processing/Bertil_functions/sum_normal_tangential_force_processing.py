import sys
import os
import re
import pandas as pd
import numpy as np


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
    normal_force_vec = []
    tangential_force_vec = []

    for i in range(0,len(time)):
        key = str(time[i])
        if time[i].is_integer():
            key = str(int(time[i]))
        ## Read file here
        file_to_open = argument_string+'/'+contact_time_and_file_name_dict[key]
        data = pd.read_csv(file_to_open).to_numpy()

        normal_force_vec.append(np.sum((data[:,6]**2)**.5))
        tangential_force_vec.append(np.sum((data[:,10]**2+data[:,11]**2+data[:,12]**2)**.5))
        time_vec.append(float(key))
    print(time_vec,normal_force_vec,tangential_force_vec)
