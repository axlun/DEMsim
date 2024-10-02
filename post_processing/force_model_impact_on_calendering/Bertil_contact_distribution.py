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


    #==CALENDERING==========================================================================================================
    # simulation_directory = '/scratch/users/axlun/DEMsim/results/article_2/final_runs_2/SN_301/2/' \
    #                    'electrode_calendering_el_pl_binder_el_pl_particle'

    #==MECHANICAL LOADING===================================================================================================
    # simulation_directory = '/scratch/users/axlun/DEMsim/results/article_2/final_runs/ref_sim/SN_1_3/electrode_mechanical_loading_hertz_tension/'

    # ==RESTING==============================================================================================================
    # simulation_directory = '/scratch/users/axlun/DEMsim/results/article_2/viscoelastic_testing/SN_1/electrode_relaxation_el_pl_binder_el_pl_particle_compression_01/'

    # ==RELAXATION==============================================================================================================
    # simulation_directory = '/scratch/users/axlun/DEMsim/results/article_2/final_runs_2/SN_2111/1/' \
    #                    'electrode_relaxation_el_pl_binder_el_pl_particle_tension'



    # simulation_directory = '/scratch/users/axlun/DEMsim/results/article_3/charge_cycling/SN_2/swelling_electrode_calendering/'
    # simulation_directory = '/scratch/users/axlun/DEMsim/results/article_3/charge_cycling/SN_2/electrode_swelling/'
    # simulation_directory = '/scratch/users/axlun/DEMsim/results/article_3/charge_cycling/SN_2/electrode_cycling_1/'
    # simulation_directory = '/scratch/users/axlun/DEMsim/results/article_3/charge_cycling/SN_2/electrode_material_scaling/'
    # simulation_directory = '/scratch/users/axlun/DEMsim/results/article_3/charge_cycling/SN_2/electrode_swelling_material_scaling/'
    # simulation_directory = '/scratch/users/axlun/DEMsim/results/article_3/charge_cycling/SN_2/swelling_electrode_mechanical_loading_compression/'
    # simulation_directory = '/scratch/users/axlun/DEMsim/results/article_3/charge_cycling/SN_2/swelling_electrode_mechanical_loading_tension/'


    # simulation_directory = '/scratch/users/axlun/DEMsim/results/article_3/charge_cycling/SN_1_reduced_cal/swelling_electrode_calendering_1115/'
    # simulation_directory = '/scratch/users/axlun/DEMsim/results/article_3/charge_cycling/SN_1_reduced_cal/swelling_electrode_calendering_1135/'

    # simulation_directory = '/scratch/users/axlun/DEMsim/results/article_3/charge_cycling/SN_3/swelling_electrode_calendering/'
    simulation_directory = '/scratch/users/axlun/DEMsim/results/article_3/charge_cycling/SN_3/electrode_swelling_material_scaling/'

    # simulation_directory = '/scratch/users/axlun/DEMsim/results/article_3/charge_cycling/SN_1/swelling_electrode_calendering/'
    # simulation_directory = '/scratch/users/axlun/DEMsim/results/article_3/result_discrepancies/SN_1/swelling_electrode_calendering/'

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

    plt.style.use('axel_style')

    time_vec, particle_contact_vec, binder_contact_vec, binder_particle_contact_vec = contact_counter_bertil(simulation_directory)




#===FIG 1 NUMBER OF CONTACTS=============================================================================================
    fig_contact_distribution, ax_contact_distribution = plt.subplots()
    lns_binder_contacts = ax_contact_distribution.plot(time_vec, binder_contact_vec, label=r'Binder contacts')
    lns_particle_contacts = ax_contact_distribution.plot(time_vec,particle_contact_vec, label=r'Particle contacts')
    lns_binder_particle_contacts = ax_contact_distribution.plot(time_vec, binder_particle_contact_vec, label=r'Binder-Particle contacts')
    ax_contact_distribution.set_ylabel("Number of contacts [-]")
    ax_contact_distribution.set_xlabel("Time [s]")
    ax_contact_distribution.legend(loc='best')
    ax_contact_distribution.set_ylim(ymin=0)
    fig_contact_distribution.tight_layout()

    fname = fig_dir + 'number_of_contacts'
    plt.savefig(fname)

    print('Plotting results')
    plt.show()

def no_of_contact_foo():
    #simulation_directory = '/scratch/users/axlun/DEMsim/results/electrode_calendaring/SN_hertz_5000p_btr_065'
    # simulation_directory = '/scratch/users/axlun/DEMsim/results/electrode_calendaring/SN00_400p_plastic_binder'


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