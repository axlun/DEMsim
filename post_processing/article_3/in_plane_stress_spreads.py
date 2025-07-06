from Bertil_functions.Bertil_functions import *
from force_model_impact_on_calendering.Bertil_mechanical_properties import stress_and_linear_strain_finder
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import shutil
import os
import numpy as np


class Simulation_data:
    def __init__(self, sim_dir, plot_contacts=0):
        self.sim_dir = sim_dir
        self.force_data, self.surface_force_index, self.surface_position_index, \
        self.surface_position_data, self.periodic_BC_data, self.force_fabric_tensor_data, \
        self.ke = bertil_data_gatherer(sim_dir)

        self.time, self.linear_strain, self.sxx, self.syy, self.szz, self.tau_xy, self.tau_xz, self.tau_yx, \
        self.tau_yz, self.tau_zx, self.tau_zy = \
            stress_and_linear_strain_finder(self.periodic_BC_data, self.force_fabric_tensor_data,
                                            self.surface_position_data)
        self.absolute_time = self.time - self.time[0]
        self.avg_in_plane_compression = - (self.sxx + self.syy) / 2

        self.height_time, self.avg_height = bertil_layer_height(sim_dir,50)

        if plot_contacts == 1:
            self.contacts_time, self.particle_contact_vec, self.binder_contact_vec, \
            self.binder_particle_contact_vec = contact_counter_bertil(sim_dir)


def in_plane_stress_processing(sim_dir, no_sims, sim_type, exclude_time, plot_contacts=0):
    simulation_dict = {}
    for i in range(1, no_sims + 1):
        simulation_dict[i] = Simulation_data(sim_dir + str(i) + sim_type, plot_contacts)

    stress_matrix = np.zeros((no_sims, simulation_dict[1].time.size))
    height_matrix = np.zeros((no_sims, simulation_dict[1].time.size))
    particle_contacts_matrix = np.zeros((no_sims, simulation_dict[1].time.size))
    binder_contacts_matrix = np.zeros((no_sims, simulation_dict[1].time.size))

    for i in range(no_sims):
        stress_matrix[i, :] = simulation_dict[i + 1].avg_in_plane_compression
        height_matrix[i, :] = simulation_dict[i + 1].avg_height
        if plot_contacts:
            particle_contacts_matrix[i, :] = simulation_dict[i + 1].particle_contact_vec
            binder_contacts_matrix[i, :] = simulation_dict[i + 1].binder_contact_vec

    time_removal_indices = np.empty(shape=[0, 1])
    time = simulation_dict[1].absolute_time
    for i in exclude_time:
        time_removal_indices = np.append(time_removal_indices, np.nonzero((time >= i[0]) & (time <= i[-1])))
    reduced_time = np.delete(time, time_removal_indices.astype(int))

    stress_mean = np.mean(stress_matrix, axis=0)
    stress_mean = np.delete(stress_mean, time_removal_indices.astype(int))

    stress_std = np.std(stress_matrix, axis=0)
    stress_std = np.delete(stress_std, time_removal_indices.astype(int))

    height_mean = np.mean(height_matrix, axis=0)
    height_mean = np.delete(height_mean, time_removal_indices.astype(int))

    height_std = np.std(height_matrix, axis=0)
    height_std = np.delete(height_std, time_removal_indices.astype(int))

    SOC = np.linspace(0, 100, stress_mean.size)

    particle_contacts_mean = np.mean(particle_contacts_matrix, axis=0)
    particle_contacts_mean = np.delete(particle_contacts_mean, time_removal_indices.astype(int))
    particle_contacts_std = np.std(particle_contacts_matrix, axis=0)
    particle_contacts_std = np.delete(particle_contacts_std, time_removal_indices.astype(int))
    binder_contacts_mean = np.mean(binder_contacts_matrix, axis=0)
    binder_contacts_mean = np.delete(binder_contacts_mean, time_removal_indices.astype(int))
    binder_contacts_std = np.std(binder_contacts_matrix, axis=0)
    binder_contacts_std = np.delete(binder_contacts_std, time_removal_indices.astype(int))

    return time, SOC, stress_mean, stress_std, stress_matrix, \
           particle_contacts_mean, particle_contacts_std, particle_contacts_matrix, \
           binder_contacts_mean, binder_contacts_std, binder_contacts_matrix, \
           height_mean, height_std, height_matrix


class Simulation:
    def __init__(self, sim_dir, no_sims, sim_type, exclude_time, plot_contacts=0):
        self.time, self.SOC, self.in_plane_stress_mean, self.in_plane_stress_std, self.stress_matrix, \
        self.particle_contacts_mean, self.particle_contacts_std, self.particle_contacts_matrix, \
        self.binder_contacts_mean, self.binder_contacts_std, self.binder_contacts_matrix, \
        self.height_mean, self.height_std, self.height_matrix \
            = in_plane_stress_processing(sim_dir, no_sims, sim_type, exclude_time, plot_contacts)


if __name__ == '__main__':
    # ==PLOT PARAMETERS=================================================================================================
    fig_dir = 'c:/temp/figures/article_3/in_plane_stress/'
    try:
        shutil.rmtree(fig_dir)
        os.mkdir(fig_dir)
    except FileNotFoundError:
        print('Directory does not exist')
        os.mkdir(fig_dir)
    except:
        print('Directory could not be removed')
        quit()
    plt.style.use('axel_style_3')
    # ==================================================================================================================

    plot_contacts = 1
    no_sims = 4
    cycling = Simulation('/scratch/users/axlun/DEMsim/results/article_3/final_runs/',
                         no_sims, '/electrode_cycling_1/',
                         [[0, 0]], plot_contacts)

    charging = Simulation('/scratch/users/axlun/DEMsim/results/article_3/final_runs/',
                          no_sims, '/electrode_swelling_material_scaling/',
                          [[0, 2], [12.30, 14.30], [23.99, 25.99]], plot_contacts)

    swelling = Simulation('/scratch/users/axlun/DEMsim/results/article_3/final_runs/',
                          no_sims, '/electrode_swelling/',
                          [[0, 2], [12.30, 14.30], [23.99, 25.99]], 0)

    mat_scaling = Simulation('/scratch/users/axlun/DEMsim/results/article_3/final_runs/',
                             no_sims, '/electrode_material_scaling/',
                             [[0, 0]], 0)


    # =PRESSURE VS SOC==================================================================================================
    fig_pressure_SOC, ax_pressure_SOC = plt.subplots()
    ax_pressure_SOC.set_xlim(xmin=0, xmax=100)
    ax_pressure_SOC.set_ylim(ymin=0, ymax=180)
    ax_pressure_SOC.set_xticks(np.linspace(0, 100, 11))
    ax_pressure_SOC.set_xlabel('SOC [%]')
    ax_pressure_SOC.set_ylabel('In-plane compressive stress [MPa]')
    ax_pressure_SOC.plot(charging.SOC, charging.in_plane_stress_mean * 1E-6, 'k-')
    ax_pressure_SOC.fill_between(charging.SOC,
                                 (charging.in_plane_stress_mean - charging.in_plane_stress_std) * 1E-6,
                                 (charging.in_plane_stress_mean + charging.in_plane_stress_std) * 1E-6,
                                 color='C0')
    fig_pressure_SOC.tight_layout()
    fname = fig_dir + 'pressure_soc_charging'
    plt.savefig(fname)
    plt.savefig(fname+'.svg')

    # =HEIGHT VS SOC====================================================================================================
    fig_height_SOC, ax_height_SOC = plt.subplots()
    ax_height_SOC.set_xlim(xmin=0, xmax=100)
    ax_height_SOC.set_ylim(ymin=0, ymax=130)
    ax_height_SOC.set_xticks(np.linspace(0, 100, 11))
    ax_height_SOC.set_yticks(np.arange(0,131,10))
    ax_height_SOC.set_xlabel('SOC [%]')
    ax_height_SOC.set_ylabel('Active layer thickness [µm]')
    ax_height_SOC.plot(charging.SOC, charging.height_mean * 1E2, 'k-')
    ax_height_SOC.fill_between(charging.SOC,
                               (charging.height_mean - charging.height_std) * 1E2,
                               (charging.height_mean + charging.height_std) * 1E2,
                               color='C0')
    fig_pressure_SOC.tight_layout()
    fname = fig_dir + 'height_soc_charging'
    plt.savefig(fname)
    plt.savefig(fname+'.svg')

    # =HEIGHT VS TIME===================================================================================================
    fig_height_time, ax_height_time = plt.subplots()
    ax_height_time.set_xlabel('Time [s]')
    ax_height_time.set_ylabel('Active layer thickness [µm]')
    # ax_height_time.set_xlim(xmin=0)
    # ax_height_time.set_ylim(ymin=0, ymax=180)
    ax_height_time.plot(charging.time, charging.height_matrix.T * 1E2)
    fig_height_time.tight_layout()
    fname = fig_dir + 'height_time'
    plt.savefig(fname)
    plt.savefig(fname+'.svg')

    # =PRESSURES VS CHARGING/SWELLING/MAT DEG===========================================================================
    fig_pressure_split, ax_pressure_split = plt.subplots()
    ax_pressure_split.set_xlabel('SOC [%]')
    ax_pressure_split.set_ylabel('In-plane compressive stress [MPa]')
    ax_pressure_split.set_xlim(xmin=0, xmax=100)
    ax_pressure_split.set_ylim(ymin=0, ymax=180)
    ax_pressure_split.set_xticks(np.linspace(0, 100, 11))

    ax_pressure_split.plot(swelling.SOC, swelling.in_plane_stress_mean * 1E-6, 'k-')
    ax_pressure_split.fill_between(swelling.SOC,
                                   (swelling.in_plane_stress_mean - swelling.in_plane_stress_std) * 1E-6,
                                   (swelling.in_plane_stress_mean + swelling.in_plane_stress_std) * 1E-6,
                                   color='C3', label='Swelling')

    ax_pressure_split.plot(charging.SOC, charging.in_plane_stress_mean * 1E-6, 'k-')
    ax_pressure_split.fill_between(charging.SOC,
                                   (charging.in_plane_stress_mean - charging.in_plane_stress_std) * 1E-6,
                                   (charging.in_plane_stress_mean + charging.in_plane_stress_std) * 1E-6,
                                   color='C0', label='Swelling + material degradation')

    ax_pressure_split.plot(mat_scaling.SOC, mat_scaling.in_plane_stress_mean * 1E-6, 'g-')
    ax_pressure_split.fill_between(mat_scaling.SOC,
                                   (mat_scaling.in_plane_stress_mean - mat_scaling.in_plane_stress_std) * 1E-6,
                                   (mat_scaling.in_plane_stress_mean + mat_scaling.in_plane_stress_std) * 1E-6,
                                   color='C2', label='Material degradation')

    ax_pressure_split.legend(loc='best')
    fig_pressure_split.tight_layout()
    fname = fig_dir + 'pressure_soc_split'
    plt.savefig(fname)
    plt.savefig(fname+'.svg')

    # =PRESSURE VS TIME=================================================================================================
    fig_pressure_time, ax_pressure_time = plt.subplots()
    ax_pressure_time.set_xlabel('Time [s]')
    ax_pressure_time.set_ylabel('In-plane compressive stress [MPa]')
    # ax_pressure_time.set_xlim(xmin=0)
    ax_pressure_time.set_ylim(ymin=0, ymax=180)
    ax_pressure_time.plot(charging.time, charging.stress_matrix.T * 1E-6)
    ax_pressure_time.plot(swelling.time, swelling.stress_matrix.T * 1E-6)
    ax_pressure_time.plot(mat_scaling.time, mat_scaling.stress_matrix.T * 1E-6)
    fig_pressure_time.tight_layout()
    fname = fig_dir + 'pressure_time_split'
    plt.savefig(fname)
    plt.savefig(fname+'.svg')

    # =CONTACTS VS SOC FOR CHARGING=====================================================================================
    particle_normalisation = 5000
    fig_contacts_SOC, ax_contacts_SOC = plt.subplots()
    ax_contacts_SOC.set_xlabel('SOC [%]')
    ax_contacts_SOC.set_ylabel('Coordination number [-]')
    ax_contacts_SOC.set_xlim(xmin=0, xmax=100)
    ax_contacts_SOC.set_xticks(np.linspace(0, 100, 11))
    ymax= 11000
    ax_contacts_SOC.set_ylim(ymin=0, ymax=ymax)
    ax_contacts_SOC.set_yticks(np.linspace(0, ymax,12))
    # ax_contacts_SOC.set_yticklabels(np.linspace(0, 2.2,12))
    ax_contacts_SOC.set_yticklabels([0,.2,.4,.6,.8,1,1.2,1.4,1.6,1.8,2,2.2])

    ax_contacts_SOC.plot(charging.SOC, charging.binder_contacts_mean, 'k-')
    ax_contacts_SOC.fill_between(charging.SOC,
                                 (charging.binder_contacts_mean - charging.binder_contacts_std),
                                 (charging.binder_contacts_mean + charging.binder_contacts_std),
                                 color='C1', label='CBD')

    ax_contacts_SOC.plot(charging.SOC, charging.particle_contacts_mean, 'k-')
    ax_contacts_SOC.fill_between(charging.SOC,
                                 (charging.particle_contacts_mean - charging.particle_contacts_std),
                                 (charging.particle_contacts_mean + charging.particle_contacts_std),
                                 color='C0', label='Particle')

    ax_contacts_SOC.legend(loc='best')
    fig_contacts_SOC.tight_layout()
    fname = fig_dir + 'contacts_soc_charging'
    plt.savefig(fname)
    plt.savefig(fname+'.svg')

    # =CONTACTS VS TIME=================================================================================================
    fig_contacts_time, ax_contacts_time = plt.subplots()
    ax_contacts_time.set_xlabel('Time [s]')
    ax_contacts_time.set_ylabel('Number of contacts [-]')

    ax_contacts_time.plot(charging.time, charging.particle_contacts_matrix.T, label='Charging P')
    ax_contacts_time.plot(charging.time, charging.binder_contacts_matrix.T, label='Charging CBD')

    ax_contacts_time.plot(swelling.time, swelling.particle_contacts_matrix.T, label='Swelling P')
    ax_contacts_time.plot(swelling.time, swelling.binder_contacts_matrix.T, label='Swelling CBD')

    ax_contacts_time.plot(mat_scaling.time, mat_scaling.particle_contacts_matrix.T, label='Mat deg P')
    ax_contacts_time.plot(mat_scaling.time, mat_scaling.binder_contacts_matrix.T, label='Mat deg CBD')
    # ax_contacts_time.legend(loc='best')
    fig_contacts_time.tight_layout()
    fname = fig_dir + 'contacts_time'
    plt.savefig(fname)
    plt.savefig(fname+'.svg')

    # =PRESSURE VS CHARGE CYCLING=======================================================================================
    fig_pressure_cycling, ax_pressure_cycling = plt.subplots()
    ax_pressure_cycling.set_xlim(xmin=0, xmax=100)
    ax_pressure_cycling.set_ylim(ymin=-30, ymax=180)
    ax_pressure_cycling.set_xlabel('SOC [%]')
    ax_pressure_cycling.set_ylabel('In-plane compressive stress [MPa]')
    ax_pressure_cycling.plot(cycling.SOC, cycling.in_plane_stress_mean * 1E-6, 'k-')
    ax_pressure_cycling.fill_between(cycling.SOC,
                                     (cycling.in_plane_stress_mean - cycling.in_plane_stress_std) * 1E-6,
                                     (cycling.in_plane_stress_mean + cycling.in_plane_stress_std) * 1E-6,
                                     color='C0')
    tick_array = np.linspace(0, 100, 21)
    tick_label_array = [0, 100, 0, 100, 0, 100, 0, 100, 0, 100, 0, 100, 0, 100, 0, 100, 0, 100, 0, 100, 0]
    ax_pressure_cycling.set_xticks(tick_array)
    ax_pressure_cycling.set_xticklabels(tick_label_array)
    # ax_pressure_cycling.minorticks_off()

    n_cycles = 10
    end_time = cycling.time[-1]
    ax_2nd_x = ax_pressure_cycling.twiny()
    ax_2nd_x.plot([0, end_time], [0, 0],linewidth=0)
    ax_2nd_x.set_xticks(np.linspace(0, end_time, n_cycles + 1), minor=False)
    ax_2nd_x.set_xticklabels(np.arange(0, n_cycles + 1, 1))
    ax_2nd_x.set_xlim(xmin=0, xmax=end_time)
    ax_2nd_x.set_xlabel('Charge cycle [-]')
    ax_pressure_cycling.tick_params(axis='x', which='minor', bottom=False, top=False)
    ax_2nd_x.tick_params(axis='x', which='minor', bottom=False, top=False)
    # ax_2nd_x.minorticks_off()

    fig_pressure_cycling.tight_layout()
    fname = fig_dir + 'pressure_soc_cycling'
    plt.savefig(fname)
    plt.savefig(fname+'.svg')

    # =HEIGHT VS CHARGE CYCLING=======================================================================================
    fig_height_cycling, ax_height_cycling = plt.subplots()
    ax_height_cycling.set_xlim(xmin=0, xmax=100)
    ax_height_cycling.set_ylim(ymin=0, ymax=130)
    ax_height_cycling.set_xlabel('SOC [%]')
    ax_height_cycling.set_ylabel('Active layer thickness [µm]')
    ax_height_cycling.plot(cycling.SOC, cycling.height_mean * 1E2, 'k-')
    ax_height_cycling.fill_between(cycling.SOC,
                                   (cycling.height_mean - cycling.height_std) * 1E2,
                                   (cycling.height_mean + cycling.height_std) * 1E2,
                                   color='C0')
    tick_array = np.linspace(0, 100, 21)
    tick_label_array = [0, 100, 0, 100, 0, 100, 0, 100, 0, 100, 0, 100, 0, 100, 0, 100, 0, 100, 0, 100, 0]
    ax_height_cycling.set_xticks(tick_array)
    # ax_height_cycling.minorticks_off()
    ax_height_cycling.set_xticklabels(tick_label_array)

    ax_2nd_x = ax_height_cycling.twiny()
    ax_2nd_x.plot([0, end_time], [0, 0],linewidth=0)
    ax_2nd_x.set_xticks(np.linspace(0, end_time, n_cycles + 1), minor=False)
    ax_2nd_x.set_xticklabels(np.arange(0, n_cycles + 1, 1))
    ax_2nd_x.set_xlim(xmin=0, xmax=end_time)
    ax_2nd_x.set_xlabel('Charge cycle [-]')
    # ax_2nd_x.minorticks_off()
    ax_height_cycling.tick_params(axis='x', which='minor', bottom=False, top=False)
    ax_2nd_x.tick_params(axis='x', which='minor', bottom=False, top=False)

    fig_height_cycling.tight_layout()
    fname = fig_dir + 'height_soc_cycling'
    plt.savefig(fname)
    plt.savefig(fname+'.svg')

    # =PRESSURE VS TIME CYCLING=========================================================================================
    fig_pressure_time, ax_pressure_time = plt.subplots()
    ax_pressure_time.set_xlabel('Time [s]')
    ax_pressure_time.set_ylabel('In-plane compressive stress [MPa]')
    ax_pressure_time.plot(cycling.time, cycling.stress_matrix.T * 1E-6)
    fig_pressure_time.tight_layout()
    fname = fig_dir + 'pressure_time_cycling'
    plt.savefig(fname)
    plt.savefig(fname+'.svg')

    # =HEIGHT VS TIME CYCLING=========================================================================================
    fig_height_time, ax_height_time = plt.subplots()
    ax_height_time.set_xlabel('Time [s]')
    ax_height_time.set_ylabel('Active layer thickness [µm]')
    ax_height_time.plot(cycling.time, cycling.height_matrix.T * 1E2)
    fig_height_time.tight_layout()
    fname = fig_dir + 'height_time_cycling'
    plt.savefig(fname)
    plt.savefig(fname+'.svg')

    # =CONTACTS VS CHARGE CYCLING=======================================================================================
    fig_contacts_cycling, ax_contacts_cycling = plt.subplots()
    ax_contacts_cycling.set_xlim(xmin=0, xmax=100)
    ax_contacts_cycling.set_xlabel('SOC [%]')
    ax_contacts_cycling.set_ylabel('Coordination number [-]')
    ymax= 11000
    ax_contacts_cycling.set_ylim(ymin=0, ymax=ymax)
    ax_contacts_cycling.set_yticks(np.linspace(0, ymax,12))
    # ax_contacts_cycling.set_yticklabels(np.linspace(0, ymax/particle_normalisation,12))
    ax_contacts_cycling.set_yticklabels([0,.2,.4,.6,.8,1,1.2,1.4,1.6,1.8,2,2.2])

    ax_contacts_cycling.plot(cycling.SOC, cycling.binder_contacts_mean, 'k-')
    ax_contacts_cycling.fill_between(cycling.SOC,
                                     (cycling.binder_contacts_mean - cycling.binder_contacts_std),
                                     (cycling.binder_contacts_mean + cycling.binder_contacts_std),
                                     color='C1', label='CBD')

    ax_contacts_cycling.plot(cycling.SOC, cycling.particle_contacts_mean, 'k-')
    ax_contacts_cycling.fill_between(cycling.SOC,
                                     (cycling.particle_contacts_mean - cycling.particle_contacts_std),
                                     (cycling.particle_contacts_mean + cycling.particle_contacts_std),
                                     color='C0', label='Particle')

    tick_array = np.linspace(0, 100, 21)
    tick_label_array = [0, 100, 0, 100, 0, 100, 0, 100, 0, 100, 0, 100, 0, 100, 0, 100, 0, 100, 0, 100, 0]
    ax_contacts_cycling.set_xticks(tick_array)
    # ax_contacts_cycling.minorticks_off()
    ax_contacts_cycling.set_xticklabels(tick_label_array)
    ax_contacts_cycling.legend(loc='best')

    ax_2nd_x = ax_contacts_cycling.twiny()
    ax_2nd_x.plot([0, end_time], [0, 0],linewidth=0)
    ax_2nd_x.set_xticks(np.linspace(0, end_time, n_cycles + 1), minor=False)
    ax_2nd_x.set_xticklabels(np.arange(0, n_cycles + 1, 1))
    ax_2nd_x.set_xlim(xmin=0, xmax=end_time)
    ax_2nd_x.set_xlabel('Charge cycle [-]')
    # ax_2nd_x.minorticks_off()
    ax_contacts_cycling.tick_params(axis='x', which='minor', bottom=False, top=False)
    ax_2nd_x.tick_params(axis='x', which='minor', bottom=False, top=False)

    fig_contacts_cycling.tight_layout()
    fname = fig_dir + 'contacts_soc_cycling'
    plt.savefig(fname)
    plt.savefig(fname+'.svg')

    # =CONTACTS VS TIME=================================================================================================
    fig_contacts_time, ax_contacts_time = plt.subplots()
    ax_contacts_time.set_xlabel('Time [s]')
    ax_contacts_time.set_ylabel('Number of contacts [-]')

    ax_contacts_time.plot(cycling.time, cycling.particle_contacts_matrix.T)
    ax_contacts_time.plot(cycling.time, cycling.binder_contacts_matrix.T)

    fig_contacts_time.tight_layout()
    fname = fig_dir + 'contacts_time_cycling'
    plt.savefig(fname)
    plt.savefig(fname+'.svg')


    print('Show plots')
    plt.show()
