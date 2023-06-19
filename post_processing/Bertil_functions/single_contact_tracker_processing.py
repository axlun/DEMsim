import sys
import os
import re
import pandas as pd
import numpy as np
np.set_printoptions(threshold=np.inf)

if __name__ == '__main__':
    argument_string = sys.argv[1]
    p1 = int(sys.argv[2])
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
    result_data_vec = []
    # contact_data_vec = []

    for i in range(0,len(time)):
        key = str(time[i])
        if time[i].is_integer():
            key = str(int(time[i]))
        ## Read file here
        file_to_open = argument_string+'/'+contact_time_and_file_name_dict[key]

        data = pd.read_csv(file_to_open).to_numpy()
        if i == 0:
            # contact_data_vec =np.zeros(len(data[0,:]))
            result_data_vec = np.array([0, 0])
            time_vec.append(float(key))
            continue
        if len(data[np.where((data[:, 0] == p1))]) != 0:
            # contact_data_vec = np.vstack([contact_data_vec,data[np.where((data[:, 0] == p1))]])

            print(type(data[np.where((data[:, 0] == p1))]))
            print(data[np.where((data[:, 0] == p1))].shape)
            print(data[np.where((data[:, 0] == p1))])
            print(type(data[np.where((data[:, 0] == p1))][0]))
            print(data[np.where((data[:, 0] == p1))][0])
            print(type(data[np.where((data[:, 0] == p1))][:,19]))
            print(data[np.where((data[:, 0] == p1))][:,19])

            # print(type(data[np.where((data[:, 0] == p1))]))
            # print(data[np.where((data[:, 0] == p1))])
            # print(type(data[np.where((data[:, 0] == p1))][:][19]))
            # print(data[np.where((data[:, 0] == p1))][:][19])
            temp_result_data_vec = np.array([np.shape((np.where((data[:, 0] == p1))))[0], np.sum(data[np.where((data[:, 0] == p1))][:,19])])
            result_data_vec = np.vstack([result_data_vec,temp_result_data_vec])
        elif len(data[np.where((data[:, 1] == p1))]) != 0:
            # contact_data_vec = np.vstack([contact_data_vec, data[np.where((data[:, 1] == p1))]])
            temp_result_data_vec = np.array([np.shape((np.where((data[:, 1] == p1))))[0], np.sum(np.where((data[:, 1] == p1))[19])])
            result_data_vec = np.vstack([result_data_vec,temp_result_data_vec])
        else:
            # contact_data_vec = np.vstack([contact_data_vec,np.zeros(len(data[0,:]))])
            result_data_vec = np.vstack([result_data_vec,np.array([0,0])])
        time_vec.append(float(key))
    print(time_vec)

    for i in range(0,len(result_data_vec[:,0])):
        for j in range(0,len(result_data_vec[0,:])):
            print(result_data_vec[i,j])

    # print(contact_data_vec )