from Bertil_functions.Bertil_functions import *
from Bertil_calendering_pressure import local_data_gatherer, bertil_data_gatherer


import matplotlib.pyplot as plt
import os
import shutil


if __name__ == '__main__':
    # ==PLOT PARAMETERS=================================================================================================
    fig_dir = 'c:/temp/figures/Bertil_layer_height/'
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
    # ==================================================================================================================

    # simulation_directory = '/scratch/users/axlun/DEMsim/results/article_3/charge_cycling/SN_1/electrode_swelling_material_scaling'
    # simulation_directory = '/scratch/users/axlun/DEMsim/results/article_3/charge_cycling/SN_2/electrode_swelling_material_scaling'
    # simulation_directory = '/scratch/users/axlun/DEMsim/results/article_3/charge_cycling/SN_1/electrode_swelling'
    # simulation_directory = '/scratch/users/axlun/DEMsim/results/article_3/charge_cycling/SN_1/electrode_material_scaling/'



    # simulation_directory = '/scratch/users/axlun/DEMsim/results/article_3/charge_cycling/SN_1/swelling_electrode_calendering/'
    # simulation_directory = '/scratch/users/axlun/DEMsim/results/article_3/result_discrepancies/SN_1/swelling_electrode_calendering/'

    # simulation_directory = '/scratch/users/axlun/DEMsim/results/article_3/charge_cycling/SN_2/swelling_electrode_calendering/'
    # simulation_directory = '/scratch/users/axlun/DEMsim/results/article_3/charge_cycling/SN_2/electrode_swelling_material_scaling'
    # simulation_directory = '/scratch/users/axlun/DEMsim/results/article_3/charge_cycling/SN_2/electrode_cycling_1/'

    simulation_directory = '/scratch/users/axlun/DEMsim/results/article_3/final_runs/2/swelling_electrode_calendering/'
    no_of_particles_in_avg = 50

    time_vec, avg_height = bertil_layer_height(simulation_directory, no_of_particles_in_avg)

    figure_avg, ax_avg = plt.subplots()
    lns_avg = ax_avg.plot(time_vec, avg_height)
    ax_avg.set_ylabel('Height [m]')
    ax_avg.set_xlabel('time [s]')
    fname = fig_dir + 'layer_height'
    plt.savefig(fname)

    figure_avg_norm, ax_avg_norm = plt.subplots()
    lns_avg = ax_avg_norm.plot(time_vec, avg_height/ avg_height[0])
    ax_avg_norm.set_ylabel('Normalised height [m/m]')
    ax_avg_norm.set_xlabel('time [s]')
    fname = fig_dir + 'normalised_layer_height'
    plt.savefig(fname)
    plt.show()


