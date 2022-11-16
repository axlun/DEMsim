from Bertil_functions.Bertil_functions import *
from Bertil_calendering_pressure import local_data_gatherer, bertil_data_gatherer

import numpy as np

import matplotlib.pyplot as plt
import re
import matplotlib
import os
import shutil
# matplotlib.style.use('classic')



if __name__ == '__main__':

#==NATUAL PACKING ======================================================================================================
#     simulation_directory = '/scratch/users/axlun/DEMsim/results/electrode_natural_packing_hertz/SN_hertz_500p_btr_8_brr_08_dt_1e0_MS_1e0_perBC_bugfix1/'
#     simulation_directory = '/scratch/users/axlun/DEMsim/results/electrode_natural_packing_hertz/SN_hertz_5000p_btr_8_brr_08_dt_1e0_MS_1e0/'
#     simulation_directory = '/scratch/users/axlun/DEMsim/results/electrode_natural_packing_hertz/SN_hertz_5000p_btr_5_brr_05_dt_1e0_MS_1e0'

#==CALENDERING==========================================================================================================
     # simulation_directory = '/scratch/users/axlun/DEMsim/results/electrode_calendering_hertz/SN_hertz_5000p_btr_8_brr_08_comp_time_20_hal_105_dt_1e2_MS_1e4'
    # simulation_directory = '/scratch/users/axlun/DEMsim/results/electrode_calendering_hertz/SN_hertz_5000p_btr_5_brr_05_comp_time_20_hal_105_dt_1e2_MS_1e4'
    simulation_directory = '/scratch/users/axlun/DEMsim/results/electrode_calendering_hertz/SN_ref_run_1_5000p_btr_5_brr_8_comp_time_20_hal_105_dt_1e2_MS_1e4'

#==MECHANICAL LOADING===================================================================================================
    # simulation_directory = '/scratch/users/axlun/DEMsim/results/electrode_mechanical_loading_hertz/SN_hertz_5000p_btr_8_brr_08_dt_1e1_MS_1e2_SR_2e-3_compression'
    # simulation_directory = '/scratch/users/axlun/DEMsim/results/electrode_mechanical_loading_hertz/SN_hertz_5000p_btr_8_brr_08_dt_1e2_MS_1e4_SR_2e-3_tension'
    # simulation_directory = '/scratch/users/axlun/DEMsim/results/electrode_mechanical_loading_hertz/SN_hertz_5000p_btr_5_brr_05_dt_5e1_MS_1e2_SR_2e-3_tension'
    #simulation_directory = '/scratch/users/axlun/DEMsim/results/electrode_mechanical_loading_hertz/SN_hertz_5000p_btr_5_brr_05_dt_5e1_MS_1e2_SR_2e-3_compression'


#==RESTING==============================================================================================================
#    simulation_directory = '/scratch/users/axlun/DEMsim/results/electrode_resting_hertz/SN_hertz_5000p_btr_8_brr_08_dt_5e1_MS_1e4_RT_10'

    time_vec, particle_contact_vec, binder_contact_vec = contact_counter_bertil(simulation_directory)

    fig_dir = 'C:/temp/figures/Bertil_contact_distribution/'
    try:
        shutil.rmtree(fig_dir)
        os.mkdir(fig_dir)
    except FileNotFoundError:
        print('Directory does not exist')
        os.mkdir(fig_dir)
    except:
        print('Directory could not be removed')
        quit()

    plt.rcParams['figure.figsize'] = (12,9)
    plt.rcParams['lines.linewidth'] = 2
    plt.rcParams['axes.labelweight'] = 'bold'
    plt.rcParams['axes.titleweight'] = 'bold'
    plt.rcParams['font.weight'] = 'bold'
    plt.rcParams['font.size'] = 20
    plt.rcParams["font.family"] = "Times New Roman"
    plt.rcParams['mathtext.default'] = 'regular'




#===FIG 1 NUMBER OF CONTACTS=============================================================================================
    fig_contact_distribution, ax_contact_distribution = plt.subplots()
    lns_binder_contacts = ax_contact_distribution.plot(time_vec, binder_contact_vec, label=r'Binder contacts')
    lnd_particle_contacts = ax_contact_distribution.plot(time_vec,particle_contact_vec, label=r'Particle contacts')
    ax_contact_distribution.set_ylabel("Number of contacts [-]")
    ax_contact_distribution.set_xlabel("Time [s]")
    ax_contact_distribution.legend(loc='best')
    ax_contact_distribution.set_ylim(ymin=0)
    fig_contact_distribution.tight_layout()

    fname = fig_dir + 'number_of_contacts'
    plt.savefig(fname)



    print('Plotting')
    plt.show()

def no_of_contact_foo():
    #simulation_directory = '/scratch/users/axlun/DEMsim/results/electrode_calendaring/SN_hertz_5000p_btr_065'
    simulation_directory = '/scratch/users/axlun/DEMsim/results/electrode_calendaring/SN00_400p_plastic_binder'


    bertil_file_output_call = 'cd '+ simulation_directory + '/contacts/' + '\nls'
    contact_files_from_bertil= input_output_return(bertil_file_output_call)
    time = []
    contact_time_and_file_name_dict  = {}
    for i in range(0,len(contact_files_from_bertil)):
        contact_file = contact_files_from_bertil[i]
        time_stamp = re.split(r'\Acontacts_',contact_file)
        time_stamp = re.split(r'.dou\n\Z', time_stamp[1])
        contact_file = re.split(r'\n\Z',contact_file)
        time.append(float(time_stamp[0]))
        contact_time_and_file_name_dict[time_stamp[0]] = contact_file[0]

    print(contact_time_and_file_name_dict)
    time.sort()
    print(time)


    for i in range(0,len(time)):
        binder_contact = 0
        particle_contact = 0
        key = str(time[i])
        if time[i].is_integer():
            key = str(int(time[i]))
        # print(key)
        contact_file_data = input_output_return('cd '+ simulation_directory + '/contacts/'+'\ncat '+contact_time_and_file_name_dict[key])
#        print(contact_file_data)
        for j in range(0,len(contact_file_data)):
            line_data = contact_file_data[j].split(', ')
            if float(line_data[7]) != 0 and float(line_data[8]) == 0.0:
                #print('Binder flag')
                binder_contact += 1
            if float(line_data[8]) != 0:
                #print('Particle flag')
                particle_contact += 1

            # if float(line_data[6]) > 0:
            #     print(line_data[6])

        print('Time: '+key)
        print('Binder contacts: '+str(binder_contact))
        print('Particle contacts: '+str(particle_contact))
        return