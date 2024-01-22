from Bertil_functions.Bertil_functions import *
from force_model_impact_on_calendering.Bertil_mechanical_properties import stress_and_linear_strain_finder
import matplotlib.pyplot as plt
import shutil
import os
# from force_model_impact_on_calendering.Bertil_calendering_pressure_multiple_simulations \
#     import calendering_plot_processing
import pandas

def relaxation_processing(Direction):

    start_index = np.argmax(np.absolute(Direction.linear_strain))
    # Normalise with the change in stress from initial stress state (compressive residual stresses)
    normalised_in_plane_stress = np.absolute((Direction.sxx[start_index:] - Direction.sxx[0]) /
                                                  (Direction.sxx[start_index]-Direction.sxx[0]))
    time = Direction.time[start_index:] - Direction.time[start_index]
    return time, normalised_in_plane_stress


class RelaxationProcessing:
    def __init__(self, direction):
        self.time, self.normalised_in_plane_stress = relaxation_processing(direction)


class Direction:
    def __init__(self, sim_dir):
        # self.time, self.calendering_surface_pressure, self.bottom_surface_pressure, self.calendering_surface_position, \
        # self.sxx, self.syy, self.szz, self.tau_xy, self.tau_xz, self.tau_yx, self.tau_yz, self.tau_zx, self.tau_zy, \
        # self.ke = calendering_plot_processing(sim_dir)

        self.force_data_compression, self.surface_force_index_compression, self.surface_position_index_compression,\
            self.surface_position_data, self.periodic_BC_data, self.force_fabric_tensor_data, \
            self.ke = bertil_data_gatherer(sim_dir)

        self.time, self.linear_strain, self.sxx, self.syy, self.szz, self.tau_xy, self.tau_xz, self.tau_yx, \
            self.tau_yz, self.tau_zx, self.tau_zy = \
            stress_and_linear_strain_finder(self.periodic_BC_data, self.force_fabric_tensor_data,
                                            self.surface_position_data)



class RelaxationSim:
    def __init__(self, sim_dir, label='Template'):
        self.sim_dir = sim_dir
        self.compression = Direction(sim_dir + '_compression')
        self.tension = Direction(sim_dir + '_tension')
        self.label = label
        self.compression_processed = RelaxationProcessing(self.compression)
        self.tension_processed = RelaxationProcessing(self.tension)


if __name__ == '__main__':
    # simulation_directory = '/scratch/users/axlun/DEMsim/results/article_2/viscoelastic_testing/SN_2_p_5t/' \
    #                        'electrode_relaxation_el_pl_binder_el_pl_particle_eps_dot_2e_2'
    # simulation_directory = '/scratch/users/axlun/DEMsim/results/article_2/final_runs_2/SN_2112/1/' \
    #                        'electrode_relaxation_el_pl_binder_el_pl_particle'
    # simulation_directory = '/scratch/users/axlun/DEMsim/results/article_2/final_runs_2/SN_2011/1/' \
    #                         'electrode_relaxation_el_pl_binder_el_pl_particle'
    # simulation_directory = '/scratch/users/axlun/DEMsim/results/article_2/final_runs_2/SN_2011/2/' \
    #                        'electrode_relaxation_el_pl_binder_el_pl_particle'
    # simulation_directory = '/scratch/users/axlun/DEMsim/results/article_2/final_runs_2/SN_2011/3/' \
    #                        'electrode_relaxation_el_pl_binder_el_pl_particle'
    # simulation_directory = '/scratch/users/axlun/DEMsim/results/article_2/final_runs_2/SN_2011/4/' \
    #                        'electrode_relaxation_el_pl_binder_el_pl_particle'


    simulation_directory = '/scratch/users/axlun/DEMsim/results/article_2/final_runs_2/SN_2111/1/' \
                           'electrode_relaxation_el_pl_binder_el_pl_particle'
    # simulation_directory = '/scratch/users/axlun/DEMsim/results/article_2/final_runs_2/SN_2111/2/' \
    #                        'electrode_relaxation_el_pl_binder_el_pl_particle'
    # simulation_directory = '/scratch/users/axlun/DEMsim/results/article_2/final_runs_2/SN_2111/3/' \
    #                        'electrode_relaxation_el_pl_binder_el_pl_particle'
    # simulation_directory = '/scratch/users/axlun/DEMsim/results/article_2/final_runs_2/SN_2111/4/' \
    #                        'electrode_relaxation_el_pl_binder_el_pl_particle'

    # simulation_directory = '/scratch/users/axlun/DEMsim/results/article_2/final_runs_2/SN_2112/3/' \
    #                        'electrode_relaxation_el_pl_binder_el_pl_particle'

    simulation = RelaxationSim(simulation_directory)

    # =EXPERIMENTAL RESULTS=============================================================================================
    #df_experimental_tension = pandas.read_excel(
    #   io='G:\My Drive\Skola\KTH\PhD\Litteratur\Characterization of the Constitutive Behavior of a Cathode Active Layer '
    #      'in Lithium-Ion Batteries Using a Bending Test Method/ExperimentData_RefinedNew_ToAxel.xlsx',
    #   sheet_name='Tension15',
    #   index_col=None)
    #experimental_tension_time = df_experimental_tension['Time'].to_numpy()[1:]
    #experimental_tension_force = df_experimental_tension['Standard force (M3)'].to_numpy()[1:]

    df_experimental_tension = pandas.read_csv(
        filepath_or_buffer='G:\My Drive\Skola\KTH\PhD\Litteratur\Characterization of the Constitutive Behavior of a '
                           'Cathode Active Layer in Lithium-Ion Batteries Using a Bending Test Method/'
                           'Relaxation_tension.csv')
    experimental_tension_time = df_experimental_tension['Time [s]'].to_numpy()
    experimental_tension_force = df_experimental_tension['Force [N]'].to_numpy()

    #df_experimental_compression= pandas.read_excel(
    #    io='G:\My Drive\Skola\KTH\PhD\Litteratur\Characterization of the Constitutive Behavior of a Cathode Active Layer '
    #       'in Lithium-Ion Batteries Using a Bending Test Method/ExperimentData_RefinedNew_ToAxel.xlsx',
    #    sheet_name='Compression15',
    #    index_col=None)
    #experimental_compression_time = df_experimental_tension['Time'].to_numpy()[1:]
    #experimental_compression_force = df_experimental_tension['Standard force (M3)'].to_numpy()[1:]

    df_experimental_compression= pandas.read_csv(
        filepath_or_buffer='G:\My Drive\Skola\KTH\PhD\Litteratur\Characterization of the Constitutive Behavior of a '
                           'Cathode Active Layer in Lithium-Ion Batteries Using a Bending Test Method/'
                           'Relaxation_compression.csv')
    experimental_compression_time = df_experimental_compression['Time [s]'].to_numpy()
    experimental_compression_force = df_experimental_compression['Force [N]'].to_numpy()

    # =PLOTTING PARAMETERS==============================================================================================
    fig_dir = 'c:/temp/figures/Bertil_relaxation/'

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
    # ==FIGURE 1 ALL STRESS TENSION=====================================================================================
    fig_all_stress_time_tension, ax_all_stress_time_tension = plt.subplots()
    lns_sxx = ax_all_stress_time_tension.plot(simulation.tension.time,
                                              simulation.tension.sxx * 1E-6, label=r'$\sigma_{xx}$')
    lns_syy = ax_all_stress_time_tension.plot(simulation.tension.time,
                                              simulation.tension.syy * 1E-6, label=r'$\sigma_{yy}$')
    lns_szz = ax_all_stress_time_tension.plot(simulation.tension.time,
                                              simulation.tension.szz * 1E-6, label=r'$\sigma_{zz}$')
    lns_txy = ax_all_stress_time_tension.plot(simulation.tension.time,
                                              simulation.tension.tau_xy * 1E-6, label=r'$\tau_{xy}$')
    lns_txz = ax_all_stress_time_tension.plot(simulation.tension.time,
                                              simulation.tension.tau_xz * 1E-6, label=r'$\tau_{xz}$')
    lns_tyx = ax_all_stress_time_tension.plot(simulation.tension.time,
                                              simulation.tension.tau_yx * 1E-6, label=r'$\tau_{yx}$')
    lns_tzx = ax_all_stress_time_tension.plot(simulation.tension.time,
                                              simulation.tension.tau_zx * 1E-6, label=r'$\tau_{zx}$')
    lns_tzy = ax_all_stress_time_tension.plot(simulation.tension.time,
                                              simulation.tension.tau_zy * 1E-6, label=r'$\tau_{zy}$')

    ax_all_stress_time_tension_2 = ax_all_stress_time_tension.twinx()
    lns_tension = ax_all_stress_time_tension_2.plot(simulation.tension.time,
                                                    simulation.tension.linear_strain * 1E2, '--')
    ax_all_stress_time_tension_2.set_ylabel('Strain [%]')
    ax_all_stress_time_tension.legend(loc='best')
    ax_all_stress_time_tension.set_xlabel('Time [s]')
    ax_all_stress_time_tension.set_ylabel('Stress [MPa]')
    ax_all_stress_time_tension.set_title('Relaxation in tension')
    fname = fig_dir + 'stress_time_tension'
    plt.savefig(fname)

    fig_ke_tension, ax_ke_tension = plt.subplots()
    ax_ke_tension.plot(simulation.tension.time, simulation.tension.ke[:, 2])
    ax_ke_tension.set_xlabel('Time [s]')
    ax_ke_tension.set_ylabel('Kinetic energy [J]')
    ax_ke_tension.set_title('Kinetic energy in tension')
    fname = fig_dir + 'ke_tension'
    plt.savefig(fname)

    # ==FIGURE 2 ALL STRESS COMPRESSION=================================================================================
    fig_all_stress_time_compression, ax_all_stress_time_compression = plt.subplots()
    lns_sxx = ax_all_stress_time_compression.plot(simulation.compression.time,
                                                  simulation.compression.sxx * 1E-6, label=r'$\sigma_{xx}$')
    lns_syy = ax_all_stress_time_compression.plot(simulation.compression.time,
                                                  simulation.compression.syy * 1E-6, label=r'$\sigma_{yy}$')
    lns_szz = ax_all_stress_time_compression.plot(simulation.compression.time,
                                                  simulation.compression.szz * 1E-6, label=r'$\sigma_{zz}$')
    lns_txy = ax_all_stress_time_compression.plot(simulation.compression.time,
                                                  simulation.compression.tau_xy * 1E-6, label=r'$\tau_{xy}$')
    lns_txz = ax_all_stress_time_compression.plot(simulation.compression.time,
                                                  simulation.compression.tau_xz * 1E-6, label=r'$\tau_{xz}$')
    lns_tyx = ax_all_stress_time_compression.plot(simulation.compression.time,
                                                  simulation.compression.tau_yx * 1E-6, label=r'$\tau_{yx}$')
    lns_tzx = ax_all_stress_time_compression.plot(simulation.compression.time,
                                                  simulation.compression.tau_zx * 1E-6, label=r'$\tau_{zx}$')
    lns_tzy = ax_all_stress_time_compression.plot(simulation.compression.time,
                                                  simulation.compression.tau_zy * 1E-6, label=r'$\tau_{zy}$')
    ax_all_stress_time_compression_2 = ax_all_stress_time_compression.twinx()
    lns_tension = ax_all_stress_time_compression_2.plot(simulation.compression.time,
                                                        simulation.compression.linear_strain * 1E2, '--')
    ax_all_stress_time_compression_2.set_ylabel('Strain [%]')
    ax_all_stress_time_compression.legend(loc='best')
    ax_all_stress_time_compression.set_xlabel('Time [s]')
    ax_all_stress_time_compression.set_ylabel('Stress [MPa]')
    ax_all_stress_time_compression.set_title('Relaxation in compression')
    fname = fig_dir + 'stress_time_compression'
    plt.savefig(fname)

    fig_ke_compression, ax_ke_compression = plt.subplots()
    ax_ke_compression.plot(simulation.compression.time, simulation.compression.ke[:, 2])
    ax_ke_compression.set_xlabel('Time [s]')
    ax_ke_compression.set_ylabel('Kinetic energy [J]')
    ax_ke_compression.set_title('Kinetic energy in compression')
    fname = fig_dir + 'ke_compression'
    plt.savefig(fname)

    # =NORMALISED DATA TENSION==========================================================================================
    fig_normalised_tension, ax_normalised_tension = plt.subplots()
    lns_experiment_normalised = ax_normalised_tension.plot(
        experimental_tension_time, experimental_tension_force/experimental_tension_force[0],
        label='Experiment')
    lns_normalised_sxx = ax_normalised_tension.plot(simulation.tension_processed.time*1E3,
                                                    simulation.tension_processed.normalised_in_plane_stress,
                                                    label=simulation.label)
    ax_normalised_tension.set_ylabel('Normalised stress [MPa/MPa]')
    ax_normalised_tension.set_xlabel('Time [s]')
    ax_normalised_tension.set_title('Relaxation in tension')
    ax_normalised_tension.set_ylim(ymin=0)
    ax_normalised_tension.set_xlim(xmin=-60)
    ax_normalised_tension.legend(loc='best')
    fname = fig_dir + 'normalised_stress_time_tension'
    plt.savefig(fname)


    # =NORMALISED DATA COMPRESSION======================================================================================
    fig_normalised_compression, ax_normalised_compression = plt.subplots()
    lns_experiment_normalised = ax_normalised_compression.plot(
        experimental_compression_time, experimental_compression_force/experimental_compression_force[0],
        label='Experiment')
    lns_normalised_sxx = ax_normalised_compression.plot(simulation.compression_processed.time*1E3,
                                                        simulation.compression_processed.normalised_in_plane_stress,
                                                        label=simulation.label)
    ax_normalised_compression.set_ylabel('Normalised stress [MPa/MPa]')
    ax_normalised_compression.set_xlabel('Time [s]')
    ax_normalised_compression.set_title('Relaxation in compression')
    ax_normalised_compression.set_ylim(ymin=0)
    ax_normalised_compression.set_xlim(xmin=-60)
    ax_normalised_compression.legend(loc='best')
    fname = fig_dir + 'normalised_stress_time_compression'
    plt.savefig(fname)

    plt.show()
