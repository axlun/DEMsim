
import paramiko
import sys
sys.path.insert(0, 'c:/Users/Axel/Documents/Secrets')
from secrets import pw

import numpy


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


def commad_input(command0):
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
    return

if __name__ == '__main__':
    command0 = "cd /scratch/users/axlun/DEMsim/"
    command = "ls"

    input_command = command0 +" \n " + command
    commad_input(input_command)
