from Bertil_functions.Bertil_functions import *

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from os.path import exists
import shutil
import os
from io import StringIO
import pandas as pd



if __name__ == '__main__':
    plt.style.use('axel_style')

    # =PACKING==========================================================================================================
    simulation_directory = '/scratch/users/axlun/DEMsim/results/article_2/final_runs/particle_contact_model/SN_2/electrode_natural_packing_el_pl_binder_el_pl_particle/'
    time_step = 12.12

    # =CALENDERING======================================================================================================
    # simulation_directory = '/scratch/users/axlun/DEMsim/results/final_runs/SN_run_1/electrode_calendering_hertz'
    # time_step = 32.5975

    # simulation_directory = '/scratch/users/axlun/DEMsim/results/article_2/contact_model_testing/SN_3/electrode_calendering_el_pl_binder_el_pl_particle/'
    # time_step = 32.5981

    # simulation_directory = '/scratch/users/axlun/DEMsim/results/final_runs/SN_run_1_El_Pl/electrode_calendering_hertz'
    # time_step = 32.5978

    # ==================================================================================================================

    # =MECHANICAL LOADING===============================================================================================
    # simulation_directory = '/scratch/users/axlun/DEMsim/results/final_runs/SN_run_1/electrode_mechanical_loading_hertz_compression'
    # time_step = 59.3712

    # ==================================================================================================================
    if (one_file_reader(simulation_directory+'/particles/particles_'+str(time_step)+'.dou') == OSError):
        print("Particle file not existing on Bertil")
        exit()

    particle_string = ''.join(one_file_reader(simulation_directory + '/particles/particles_' + str(time_step) + '.dou'))
    particle_file = StringIO(particle_string)
    particle_df = pd.read_csv(particle_file,
                              names=['Id','x','y','z','Rot_x','Rot_y','Rot_z','R','E_k_p','Mat-id','F_x', 'F_y','F_z'])

    overlap_measure = 'h_'
    # print(particle_df)

    contact_string = ''.join(one_file_reader(simulation_directory+'/contacts/contacts_'+str(time_step)+'.dou'))
    contact_file = StringIO(contact_string)
    contact_df = pd.read_csv(contact_file,
                             names=['P_1','P_2','n_x','n_y','n_z','h_','F_','F_b','F_p','h_max','F_Tx','F_Ty','F_Tz',
                                    'F_Tbx','F_Tby','F_Tbz','F_Tpx','F_Tpy','F_Tpz','Binder contact','b_t'])

    particle_contact_df = contact_df[contact_df[overlap_measure] > 0]
    particle_contact_df = particle_contact_df[particle_contact_df['P_1'] > particle_df.iloc[0]['Id']]
    particle_contact_df = particle_contact_df[particle_contact_df['P_2'] > particle_df.iloc[0]['Id']]

    hist_vec = []
    for index, row in particle_contact_df.iterrows():
        particle_1 = particle_df[particle_df['Id'] == row['P_1']]['R']
        particle_1 = particle_1.to_numpy()[0]
        particle_2 = particle_df[particle_df['Id'] == row['P_2']]['R']
        particle_2 = particle_2.to_numpy()[0]
        R_eff = particle_1*particle_2/(particle_1+particle_2)
        h_ = row[overlap_measure]
        hist_vec.append(h_/R_eff)
        if h_/R_eff > 1:
            print(row['P_1'],particle_1)
            print(row['P_2'],particle_2)
            print(h_)

    # print(hist_vec)
    yield_point = 0.00859
    yield_point = 0.0729

    bins = np.arange(start=0, stop=2, step=yield_point/2)
    hist, bin_edges=np.histogram(hist_vec,bins=bins)
    print(np.sum(hist))

    fig,ax = plt.subplots()
    lns = ax.hist(hist_vec,bins=bins,density=False)
    lns_yield_point = ax.plot([yield_point,yield_point],[0,1.1*max(hist)],label='Yield overlap')
    ax.set_xlim(xmin=0,xmax=1)
    ax.set_ylim(ymin=0,ymax=1.1*max(hist))
    ax.set_xlabel('$h/R^*$')
    ax.set_ylabel('Number of contacts')
    ax.legend(loc='best')
    # lns_data_hist = ax_data_hist.hist(Ebner_q_0_sample,bins=100,density=True)


    plt.show()
    # print(particle_file[-1])
    # print(contact_file[-1])
    #
    # print(particle_file[-1][:-1])
    # print(contact_file[-1][:-1])


