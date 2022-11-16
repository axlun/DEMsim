import os

import numpy as np
import re
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.style
import time
import pandas as pd




if __name__ == '__main__':

    # ==NATUAL PACKING ======================================================================================================
    # simulation_directory = "c:/Users/Axel/Documents/DEM/results/electrode_natural_packing_hertz/SN_hertz_200p_btr_8_brr_08_dt_1e0_MS_1e0_elast_binder_new_tang/"
    # simulation_directory = "c:/Users/Axel/Documents/DEM/results/electrode_natural_packing_hertz/SN_hertz_200p_btr_8_brr_08_dt_1e0_MS_1e0_elast_binder_new_tang_W_perBC/"
    simulation_directory = "c:/Users/Axel/Documents/DEM/Bertil_results/electrode_natural_packing_hertz/SN_hertz_5000p_btr_8_brr_08_dt_1e0_MS_1e0_el_b_new_tang_no_PerBC/"
    # simulation_directory = "c:/Users/Axel/Documents/DEM/Bertil_results/electrode_natural_packing_hertz/SN_hertz_5000p_btr_8_brr_08_dt_1e0_MS_1e0_el_b_new_tang/"


    #==CALENDERING======================================================================================================
    # simulation_directory =  "C:/Users/Axel/Documents/DEM/Bertil_results/electrode_calendering_hertz/SN_hertz_5000p_btr_8_brr_08_comp_time_40_hal_105_dt_1e2_MS_1e4_from_larger_particles/"

    # ==MECHANICAL LOADING==============================================================================================
    # simulation_directory =  "C:/Users/Axel/Documents/DEM/Bertil_results/electrode_mechanical_loading_hertz/SN_hertz_5000p_btr_8_brr_08_large_part_dt_5e1_MS_1e4_SR_1e-3_compression/"
    # simulation_directory =  "C:/Users/Axel/Documents/DEM/results/contact_testing/elastic_plastic_binder_hertz_particle/tangential_force/New_tangential_force_relation/"



    contact_files = os.listdir(simulation_directory+"contacts/")
    time_steps = []
    particle_time_and_file_name_dict = {}
    for i in range(0,len(contact_files)):
        contact_file = contact_files[i]
        time_stamp = re.split(r'\Acontacts_', contact_file)
        time_stamp = re.split(r'.dou\Z', time_stamp[1])
        time_steps.append(float(time_stamp[0]))
        particle_time_and_file_name_dict[time_stamp[0]] = contact_file
    time_steps.sort()

    contact_data_dict = {}
    normal_force = []
    tangential_force = []
    number_of_points = int(len(time_steps)/1)
    print(number_of_points)
    for i in range(0, number_of_points):#len(time)):
        start_time = time.time()
        key = str(time_steps[i])
        if time_steps[i].is_integer():
            key = str(int(time_steps[i]))
        data_read_time = time.time()
        data = pd.read_csv(simulation_directory + 'contacts/' + particle_time_and_file_name_dict[key]).to_numpy()
        # data = np.genfromtxt(simulation_directory + 'contacts/' + particle_time_and_file_name_dict[key],
        #                      delimiter=', ')
        # print(time_steps[i])
        # print('Data read time = '+str(time.time()-data_read_time)+'s \n')

        contact_data_dict[time_steps[i]] = data
        normal_force.append(np.sum((data[:,6]**2)**.5))
        tangential_force.append(np.sum((data[:,10]**2+data[:,11]**2+data[:,12]**2)**.5))

        # print('Execute time = '+str(time.time()-start_time)+'s \n')




    # ========NORMAL AND TANGENTIAL FORCES==============================================================================
    print("Plotting results")
    figure_normal_tangential_force_time, ax_normal_tangential_force_time = plt.subplots()
    lns_normal_force_time = ax_normal_tangential_force_time.plot(time_steps[:number_of_points], normal_force,'r',linewidth=3,label=r'F_n')
    lns_tangential_force_time = ax_normal_tangential_force_time.plot(time_steps[:number_of_points], tangential_force,'g',linewidth=3,label=r'F_t')
    ax_normal_tangential_force_time.set_xlabel('Time [s]')
    ax_normal_tangential_force_time.set_ylabel('Total force[N]')
    ax_normal_tangential_force_time.legend()
    plt.show()
