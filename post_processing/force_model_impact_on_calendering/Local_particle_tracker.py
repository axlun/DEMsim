from Bertil_functions.Bertil_functions import *
from Bertil_calendering_pressure import local_data_gatherer, bertil_data_gatherer

import numpy as np

import matplotlib.pyplot as plt
import sys
import os
import re
import matplotlib
import pandas as pd

matplotlib.style.use('classic')

def particle_tracker_local(sim_dir,p1):

    contact_files = os.listdir(sim_dir+'/contacts')
    particle_files = os.listdir(sim_dir+'/particles')
    time = []
    particle_time_and_file_name_dict = {}
    for i in range(0,len(particle_files)):
        particle_file = particle_files[i]
        time_stamp = re.split(r'\Aparticles_',particle_file)
        time_stamp = re.split(r'.dou\Z', time_stamp[1])

        time.append(float(time_stamp[0]))
        particle_time_and_file_name_dict[time_stamp[0]] = particle_file

    time.sort()
    time_vec = []
    particle_data_vec = []
    for i in range(0,len(time)):
        key = str(time[i])
        if time[i].is_integer():
            key = str(int(time[i]))
        ## Read file here
        file_to_open = sim_dir+'particles/'+particle_time_and_file_name_dict[key]

        data = pd.read_csv(file_to_open).to_numpy()

        if i == 0:
            particle_data_vec = data[np.where(data[:, 0] == p1)]
        else:
            particle_data_vec = np.vstack([particle_data_vec,data[np.where(data[:, 0] == p1)]])
        time_vec.append(float(key))
    time_vec = np.array(time_vec, dtype=float)
    return time_vec,particle_data_vec

if __name__ == '__main__':

#    simulation_directory = 'c:/Users/Axel/Documents/DEM/results/electrode_natural_packing/Build_8/'
#     simulation_directory = "c:/Users/Axel/Documents/DEM/results/electrode_natural_packing_hertz/SN_hertz_200p_btr_8_brr_08_dt_1e0_MS_1e0_elast_binder_new_tang/"

#    simulation_directory = 'c:/Users/Axel/Documents/DEM/results/natural_packing_hertz/meter_particles_N_200_mass_sclaing_1e0_gravity_1e1_time_step_1e_test1/'
    simulation_directory = 'c:/Users/Axel/Documents/DEM/Bertil_results/electrode_natural_packing_hertz/SN_hertz_2000p_btr_8_brr_08_dt_1e0_MS_1e0_el_b_new_tang_fix_gate3/'

    particle = 735

    time_vec, particle_data_vec = particle_tracker_local(simulation_directory,particle)
    print(time_vec)
    print(particle_data_vec[:,0])
    print(len(time_vec))
    print(len(particle_data_vec))


    fig_particle_force, ax_particle_force = plt.subplots()
    lns_total_force = ax_particle_force.plot(time_vec,(particle_data_vec[:,10]**2+particle_data_vec[:,11]**2+particle_data_vec[:,12]**2)**.5, 'g', label=r'Total')
    ax_particle_force.set_ylabel("Force [N]")
    ax_particle_force.set_xlabel("time [s]")
    ax_particle_force.legend(loc='best')


    fig_ke, ax_ke = plt.subplots()
    lns_h_ = ax_ke.plot(time_vec,particle_data_vec[:,8], 'r', label=r'KE')
    ax_ke.set_ylabel("Kinetic energy [J]")
    ax_ke.set_xlabel("time [s]")
    ax_ke.legend(loc='best')

    plt.show()