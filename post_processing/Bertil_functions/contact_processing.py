import sys
import os
import re


if __name__ == '__main__':
    argument_string = sys.argv[1]

    print('Hello world')

    print(argument_string)
    contact_files = os.listdir(argument_string)
#    print(contact_files)
    time = []
    contact_time_and_file_name_dict  = {}
    for i in range(0,len(contact_files)):
        contact_file = contact_files[i]
        time_stamp = re.split(r'\Acontacts_',contact_file)
        time_stamp = re.split(r'.dou\n\Z', time_stamp[1])
        contact_file = re.split(r'\n\Z',contact_file)
        time.append(float(time_stamp[0]))
        contact_time_and_file_name_dict[time_stamp[0]] = contact_file[0]

    print(contact_time_and_file_name_dict)
    time.sort()
    print(time)