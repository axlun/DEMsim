import sys
import os
import re
import pandas as pd
import numpy as np
np.set_printoptions(threshold=np.inf)

if __name__ == '__main__':
    argument_string = sys.argv[1]
    p1 = int(sys.argv[2])
    particle_files = os.listdir(argument_string)
    time = []
    particle_time_and_file_name_dict = {}
    for i in range(0, len(particle_files)):
        particle_file = particle_files[i]
        time_stamp = re.split(r'\Aparticles_', particle_file)
        time_stamp = re.split(r'.dou\Z', time_stamp[1])

        time.append(float(time_stamp[0]))
        particle_time_and_file_name_dict[time_stamp[0]] = particle_file
    time.sort()
    time_vec = []
    result_data_vec = []
    for i in range(0, len(time)):
        key = str(time[i])
        if time[i].is_integer():
            key = str(int(time[i]))
        ## Read file here
        file_to_open = argument_string + '/' + particle_time_and_file_name_dict[key]
        temp_result_data_vec = np.array([0.0, 0.0])
        data = pd.read_csv(file_to_open,names=list(range(13)))
        particle_data = data.loc[data[0] == p1].to_numpy()
        if i == 0:
            time_vec.append(float(key))
            result_data_vec = particle_data
            continue
        result_data_vec = np.vstack([result_data_vec, particle_data])
        time_vec.append(float(key))
    print(time_vec)
    print(np.array2string(result_data_vec.flatten(), separator=',', max_line_width=int(1e99)))