from Bertil_functions.Bertil_functions import *
import matplotlib.pyplot as plt
import shutil
import os
from force_model_impact_on_calendering.Bertil_calendering_pressure_multiple_simulations \
    import calendering_plot_processing
import pandas
import numpy as np

from stress_relaxation import Direction#, relaxation_processing, RelaxationProcessing


def relaxation_processing(Direction):
    start_index = np.argmax(np.absolute(Direction.linear_strain))
    # Normalise with the change in stress from initial stress state (compressive residual stresses)
    in_plane_stress = np.absolute(Direction.sxx[start_index:])# - Direction.sxx[0])
    time = Direction.time[start_index:] - Direction.time[start_index]
    return time, in_plane_stress


class RelaxationProcessing:
    def __init__(self, direction):
        self.time, self.in_plane_stress = relaxation_processing(direction)


def relaxation_spread_processing(sim_dir, sim_type, no_sims, direction):
    simulation_directory = {}
    for i in range(1,no_sims+1):
        simulation_directory[i] = RelaxationProcessing(Direction(sim_dir + str(i) + '/electrode_relaxation_' +
                                                                 sim_type + '_' + direction))
    stress_matrix = np.zeros((no_sims, simulation_directory[1].time.size))
    for i in range(no_sims):
        stress_matrix[i,:] = simulation_directory[i+1].in_plane_stress
    time = simulation_directory[1].time - simulation_directory[1].time[0]
    stress_mean = np.mean(stress_matrix, axis=0)
    stress_std = np.std(stress_matrix, axis=0)
    return time, stress_mean, stress_std, stress_matrix


class Relaxation:
    def __init__(self, sim_dir, sim_type, no_sims, direction):
        self.time, self.stress_mean, self.stress_std, self.stress_matrix = \
            relaxation_spread_processing(sim_dir, sim_type, no_sims, direction)


class Simulation:
    def __init__(self, sim_dir, sim_type, no_sims, label='Template'):
        self.sim_dir = sim_dir
        self.sim_type = sim_type
        self.no_sims = no_sims
        self.compression = Relaxation(sim_dir, sim_type, no_sims, 'compression')
        self.tension = Relaxation(sim_dir, sim_type, no_sims, 'tension')
        self.label = label

if __name__=='__main__':
    simulation_directory_SN_2011 = '/scratch/users/axlun/DEMsim/results/article_2/final_runs_2/SN_2011/'
    SN_2011 = Simulation(simulation_directory_SN_2011, 'el_pl_binder_el_pl_particle', 4, 'Prony series 1')


    simulation_directory_SN_2111 = '/scratch/users/axlun/DEMsim/results/article_2/final_runs_2/SN_2111/'
    SN_2111 = Simulation(simulation_directory_SN_2111, 'el_pl_binder_el_pl_particle', 4, 'Prony series 2')


    simulation_directory_SN_2110 = '/scratch/users/axlun/DEMsim/results/article_2/final_runs_2/SN_2110/'
    SN_2110 = Simulation(simulation_directory_SN_2110, 'el_pl_binder_el_pl_particle', 4, r'$\varepsilon = 1%$')

    simulation_directory_SN_2112 = '/scratch/users/axlun/DEMsim/results/article_2/final_runs_2/SN_2112/'
    SN_2112 = Simulation(simulation_directory_SN_2112, 'el_pl_binder_el_pl_particle', 4, r'$\varepsilon = 2%$')

    # =EXPERIMENTAL RESULTS=============================================================================================
    df_experimental_tension = pandas.read_csv(
        filepath_or_buffer='G:\My Drive\Skola\KTH\PhD\Litteratur\Characterization of the Constitutive Behavior of a '
                           'Cathode Active Layer in Lithium-Ion Batteries Using a Bending Test Method/'
                           'Relaxation_tension.csv')
    experimental_tension_time = df_experimental_tension['Time [s]'].to_numpy()
    experimental_tension_force = df_experimental_tension['Force [N]'].to_numpy()

    df_experimental_compression= pandas.read_csv(
        filepath_or_buffer='G:\My Drive\Skola\KTH\PhD\Litteratur\Characterization of the Constitutive Behavior of a '
                           'Cathode Active Layer in Lithium-Ion Batteries Using a Bending Test Method/'
                           'Relaxation_compression.csv')
    experimental_compression_time = df_experimental_compression['Time [s]'].to_numpy()
    experimental_compression_force = df_experimental_compression['Force [N]'].to_numpy()

    # =PLOTTING PARAMETERS==============================================================================================
    fig_dir = 'c:/temp/figures/Bertil_relaxation_spreads/'
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

    # =FIGURE TENSION===================================================================================================
    fig_tension, ax_tension = plt.subplots()
    lns_experiment_tension = ax_tension.plot(
        experimental_tension_time, experimental_tension_force/experimental_tension_force[0], label='Experiment')
    lns_SN_2011_tension_mean= ax_tension.plot(
        SN_2011.tension.time*1E3, SN_2011.tension.stress_mean/SN_2011.tension.stress_mean[0],
        'k-', zorder=1)
    # lns_SN_2011_tension_std= ax_tension.plot(
    #     SN_2011.tension.time*1E3, SN_2011.tension.stress_std)
    ax_tension.fill_between(
        SN_2011.tension.time*1E3,
        (SN_2011.tension.stress_mean - SN_2011.tension.stress_std) / SN_2011.tension.stress_mean[0],
        (SN_2011.tension.stress_mean + SN_2011.tension.stress_std) / SN_2011.tension.stress_mean[0],
        label=SN_2011.label, color='C1', alpha=1)

    lns_SN_2111_tension_mean = ax_tension.plot(
        SN_2111.tension.time*1E3, SN_2111.tension.stress_mean/SN_2111.tension.stress_mean[0],
        'k-', zorder=1)
    # lns_SN_2111_tension_std= ax_tension.plot(
    #     SN_2111.tension.time*1E3, SN_2111.tension.stress_std)
    ax_tension.fill_between(
        SN_2111.tension.time*1E3,
        (SN_2111.tension.stress_mean - SN_2111.tension.stress_std) / SN_2111.tension.stress_mean[0],
        (SN_2111.tension.stress_mean + SN_2111.tension.stress_std) / SN_2111.tension.stress_mean[0],
        label=SN_2111.label,color='C2', alpha=1)
    # ax_tension.set_title('Relaxation in tension')
    ax_tension.set_ylabel('Normalised stress [MPa/MPa]')
    ax_tension.set_xlabel('Time [s]')
    ax_tension.set_ylim(ymin=0, ymax=1)
    ax_tension.set_xlim(xmin=0, xmax=21600)
    ax_tension.legend(loc='best')
    fname = fig_dir + 'stress_time_tension'
    plt.savefig(fname)

    # =FIGURE TENSION ALL RUNS==========================================================================================
    fig_tension_2, ax_tension_2 = plt.subplots()
    ax_tension_2.plot(SN_2011.tension.time, SN_2011.tension.stress_matrix[0, :], label='SN_2011/1')
    ax_tension_2.plot(SN_2011.tension.time, SN_2011.tension.stress_matrix[1, :], label='SN_2011/2')
    ax_tension_2.plot(SN_2011.tension.time, SN_2011.tension.stress_matrix[2, :], label='SN_2011/3')
    ax_tension_2.plot(SN_2011.tension.time, SN_2011.tension.stress_matrix[3, :], label='SN_2011/4')
    ax_tension_2.plot(SN_2111.tension.time, SN_2111.tension.stress_matrix[0, :], label='SN_2111/1')
    ax_tension_2.plot(SN_2111.tension.time, SN_2111.tension.stress_matrix[1, :], label='SN_2111/2')
    ax_tension_2.plot(SN_2111.tension.time, SN_2111.tension.stress_matrix[2, :], label='SN_2111/3')
    ax_tension_2.plot(SN_2111.tension.time, SN_2111.tension.stress_matrix[3, :], label='SN_2111/4')
    ax_tension_2.set_title('Relaxation in tension')
    ax_tension_2.set_xlabel('Time [s]')
    ax_tension_2.set_ylabel('Stress [MPa]')
    ax_tension_2.legend(loc='best')
    fname = fig_dir + 'all_stress_time_tension'
    plt.savefig(fname)

    # =FIGURE COMPRESSION===============================================================================================
    fig_compression, ax_compression = plt.subplots()
    lns_experiment_compression = ax_compression.plot(
        experimental_compression_time, experimental_compression_force/experimental_compression_force[0],
        label='Experiment')
    lns_SN_2011_compression_mean = ax_compression.plot(
        SN_2011.compression.time*1E3, SN_2011.compression.stress_mean/SN_2011.compression.stress_mean[0],
        'k-', zorder=1)
    # lns_SN_2011_compression_std = ax_compression.plot(
    #     SN_2011.compression.time*1E3, SN_2011.compression.stress_std)
    ax_compression.fill_between(
        SN_2011.compression.time*1E3,
        (SN_2011.compression.stress_mean - SN_2011.compression.stress_std) / SN_2011.compression.stress_mean[0],
        (SN_2011.compression.stress_mean + SN_2011.compression.stress_std) / SN_2011.compression.stress_mean[0],
        label=SN_2011.label, color='C1', alpha=1)

    lns_SN_2111_compression_mean = ax_compression.plot(
        SN_2111.compression.time*1E3, SN_2111.compression.stress_mean/SN_2111.compression.stress_mean[0],
        'k-', zorder=1)
    # lns_SN_2111_compression_std = ax_compression.plot(
    #     SN_2111.compression.time*1E3, SN_2111.compression.stress_std)
    ax_compression.fill_between(
        SN_2111.compression.time*1E3,
        (SN_2111.compression.stress_mean - SN_2111.compression.stress_std) / SN_2111.compression.stress_mean[0],
        (SN_2111.compression.stress_mean + SN_2111.compression.stress_std) / SN_2111.compression.stress_mean[0],
        label=SN_2111.label,color='C2', alpha=1)
    # ax_compression.set_title('Relaxation in compression')
    ax_compression.set_ylabel('Normalised stress [MPa/MPa]')
    ax_compression.set_xlabel('Time [s]')
    ax_compression.set_ylim(ymin=0, ymax=1)
    ax_compression.set_xlim(xmin=0, xmax=21600)
    ax_compression.legend(loc='best')
    fname = fig_dir + 'stress_time_compression'
    plt.savefig(fname)

    # =FIGURE COMPRESSION ALL RUNS======================================================================================
    fig_compression_2, ax_compression_2 = plt.subplots()
    ax_compression_2.plot(SN_2011.compression.time, SN_2011.compression.stress_matrix[0, :], label='SN_2011/1')
    ax_compression_2.plot(SN_2011.compression.time, SN_2011.compression.stress_matrix[1, :], label='SN_2011/2')
    ax_compression_2.plot(SN_2011.compression.time, SN_2011.compression.stress_matrix[2, :], label='SN_2011/3')
    ax_compression_2.plot(SN_2011.compression.time, SN_2011.compression.stress_matrix[3, :], label='SN_2011/4')
    ax_compression_2.plot(SN_2111.compression.time, SN_2111.compression.stress_matrix[0, :], label='SN_2111/1')
    ax_compression_2.plot(SN_2111.compression.time, SN_2111.compression.stress_matrix[1, :], label='SN_2111/2')
    ax_compression_2.plot(SN_2111.compression.time, SN_2111.compression.stress_matrix[2, :], label='SN_2111/3')
    ax_compression_2.plot(SN_2111.compression.time, SN_2111.compression.stress_matrix[3, :], label='SN_2111/4')
    ax_compression_2.set_title('Relaxation in compression')
    ax_compression_2.set_xlabel('Time [s]')
    ax_compression_2.set_ylabel('Stress [MPa]')
    ax_compression_2.legend(loc='best')
    fname = fig_dir + 'all_stress_time_compression'
    plt.savefig(fname)

    # =STRAIN LEVELS====================================================================================================
    alpha = 1
    # =TENSION==========================================================================================================
    fig_strain_level_tension, ax_strain_level_tension = plt.subplots()
    ax_strain_level_tension.plot(SN_2110.tension.time*1E3,
                                 SN_2110.tension.stress_mean/SN_2110.tension.stress_mean[0], 'k-', zorder=1)
    ax_strain_level_tension.fill_between(
        SN_2110.tension.time*1E3,
        (SN_2110.tension.stress_mean - SN_2110.tension.stress_std)/SN_2110.tension.stress_mean[0],
        (SN_2110.tension.stress_mean + SN_2110.tension.stress_std)/SN_2110.tension.stress_mean[0],
        label=SN_2110.label, color='C3', zorder=1, alpha=alpha)

    ax_strain_level_tension.plot(SN_2111.tension.time*1E3,
                                 SN_2111.tension.stress_mean/SN_2111.tension.stress_mean[0], 'k-', zorder=2)
    ax_strain_level_tension.fill_between(
        SN_2111.tension.time*1E3,
        (SN_2111.tension.stress_mean - SN_2111.tension.stress_std)/SN_2111.tension.stress_mean[0],
        (SN_2111.tension.stress_mean + SN_2111.tension.stress_std)/SN_2111.tension.stress_mean[0],
        label=r'$\varepsilon=1.65$', color='C2', zorder=2, alpha=alpha)

    ax_strain_level_tension.plot(SN_2112.tension.time*1E3,
                                 SN_2112.tension.stress_mean/SN_2112.tension.stress_mean[0], 'k-', zorder=3)
    ax_strain_level_tension.fill_between(
        SN_2112.tension.time*1E3,
        (SN_2112.tension.stress_mean - SN_2112.tension.stress_std)/SN_2112.tension.stress_mean[0],
        (SN_2112.tension.stress_mean + SN_2112.tension.stress_std)/SN_2112.tension.stress_mean[0],
        label=SN_2112.label, color='C0', zorder=3, alpha=alpha)

    ax_strain_level_tension.set_xlim(xmin=0, xmax=21600)
    ax_strain_level_tension.set_ylim(ymin=0, ymax=1)
    ax_strain_level_tension.set_ylabel('Normalised stress [MPa/MPa]')
    ax_strain_level_tension.set_xlabel('Time [s]')
    ax_strain_level_tension.legend(loc='best')
    fname = fig_dir + 'strain_levels_tension'
    plt.savefig(fname)

    # =COMPRESSION======================================================================================================
    fig_strain_level_compression, ax_strain_level_compression = plt.subplots()
    ax_strain_level_compression.plot(SN_2110.compression.time*1E3,
                                     SN_2110.compression.stress_mean/SN_2110.compression.stress_mean[0],'k-', zorder=1)
    ax_strain_level_compression.fill_between(
        SN_2110.compression.time*1E3,
        (SN_2110.compression.stress_mean - SN_2110.compression.stress_std)/SN_2110.compression.stress_mean[0],
        (SN_2110.compression.stress_mean + SN_2110.compression.stress_std)/SN_2110.compression.stress_mean[0],
        label=SN_2110.label, color='C3', zorder=1, alpha=alpha)

    ax_strain_level_compression.plot(SN_2111.compression.time*1E3,
                                     SN_2111.compression.stress_mean/SN_2111.compression.stress_mean[0],'k-', zorder=2)
    ax_strain_level_compression.fill_between(
        SN_2111.compression.time*1E3,
        (SN_2111.compression.stress_mean - SN_2111.compression.stress_std)/SN_2111.compression.stress_mean[0],
        (SN_2111.compression.stress_mean + SN_2111.compression.stress_std)/SN_2111.compression.stress_mean[0],
        label=r'$\varepsilon=1.65$', color='C2', zorder=2, alpha=alpha)


    ax_strain_level_compression.plot(SN_2112.compression.time*1E3,
                                     SN_2112.compression.stress_mean/SN_2112.compression.stress_mean[0], 'k-', zorder=3)
    ax_strain_level_compression.fill_between(
        SN_2112.compression.time*1E3,
        (SN_2112.compression.stress_mean - SN_2112.compression.stress_std)/SN_2112.compression.stress_mean[0],
        (SN_2112.compression.stress_mean + SN_2112.compression.stress_std)/SN_2112.compression.stress_mean[0],
        label=SN_2112.label, color='C0', zorder=3, alpha=alpha)

    ax_strain_level_compression.set_xlim(xmin=0, xmax=21600)
    ax_strain_level_compression.set_ylim(ymin=0, ymax=1)
    ax_strain_level_compression.set_ylabel('Normalised stress [MPa/MPa]')
    ax_strain_level_compression.set_xlabel('Time [s]')
    ax_strain_level_compression.legend(loc='best')
    fname = fig_dir + 'strain_levels_compression'
    plt.savefig(fname)

    plt.show()


