from Bertil_functions.Bertil_functions import *
import matplotlib.pyplot as plt
import matplotlib
matplotlib.style.use('axel_style')


if __name__=='__main__':
    simulation_directory = '/scratch/users/axlun/DEMsim/results/article_4/particle_fracture/1/electrode_calendering_fracture/'
    simulation_time, number_of_fractures, fracture_array = fractured_particle_gatherer(simulation_directory)

    figure_fracture, ax_fracture = plt.subplots()
    ax_fracture.plot(simulation_time, number_of_fractures)
    ax_fracture.set_xlabel('Time [s]')
    ax_fracture.set_ylabel('Number of fractured particles [-]')
    figure_fracture.tight_layout()

    numer_of_particles = 2500

    figure_fracture, ax_fracture_normalised = plt.subplots()
    ax_fracture_normalised.plot(simulation_time, 1E2*number_of_fractures/numer_of_particles)
    ax_fracture_normalised.set_xlabel('Time [s]')
    ax_fracture_normalised.set_ylabel('Percentage of fractured particles [-]')
    ax_fracture_normalised.set_ylim(ymin=0, ymax=100)
    ax_fracture_normalised.set_yticks(range(0,101,10))
    figure_fracture.tight_layout()

    plt.show()