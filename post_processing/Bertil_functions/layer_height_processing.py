import sys
import os
import re
import pandas as pd
import numpy as np

np.set_printoptions(threshold=np.inf)

if __name__ == '__main__':
    argument_string = sys.argv[1]
    particle_files = os.listdir(argument_string)
    time = []
    particle_time_and_file_name_dict = {}
    n_particles = 0
    for i in range(0, len(particle_files)):
        particle_file = particle_files[i]
        time_stamp = re.split(r'\Aparticles_', particle_file)
        time_stamp = re.split(r'.dou\Z', time_stamp[1])
        time.append(float(time_stamp[0]))
        particle_time_and_file_name_dict[time_stamp[0]] = particle_file
    time.sort()
    time_vec = []
    avg_height_vec = np.zeros(shape=(int(len(time)) - 1))
    n_in_avg = 10
    for i in range(0, len(time)):
        key = str(time[i])
        if time[i].is_integer():
            key = str(int(time[i]))
        ## Read file here
        file_to_open = argument_string + '/' + particle_time_and_file_name_dict[key]
        # data = pd.read_csv(file_to_open, names=list(range(13)))
        data = pd.read_csv(file_to_open,header=None).to_numpy()
        avg_height_vec[i] = np.sum(np.sort(data[:,3] + data[:,7])[-n_in_avg:])/n_in_avg
        time_vec.append(float(key))
    print(time_vec)
    print(np.array2string(avg_height_vec.flatten(), separator=',', max_line_width=int(1e99)))
