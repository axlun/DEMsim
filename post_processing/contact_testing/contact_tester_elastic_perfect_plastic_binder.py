import numpy as np

import matplotlib.pyplot as plt
import matplotlib
import matplotlib.style


def get_parameter(parameter_dir, parameter):
    parameter_array = np.genfromtxt(parameter_dir, dtype='U30,f', delimiter='=')
    for index, item in np.ndenumerate(parameter_array):
        if parameter_array[index]['f0'] == parameter:
            return parameter_array[index]['f1']


def main():
    matplotlib.style.use('classic')
    plt.rc('font', serif='Computer Modern Roman')
    plt.rcParams.update({'font.size': 15})
    plt.rcParams['lines.linewidth'] = 2
    force_data = np.genfromtxt('C:/Users/Axel/Documents/DEM/results/Elastic_perfect_plastic_binder_elastic_particle_contact_test/contact_testing.dou',
                                delimiter=',')
    simulation_dir = 'C:/Users/Axel/Documents/DEM/DEMsim/simulations/Elastic_perfect_plastic_binder_contact_test/sim_1.sim'


    time_parameter = get_parameter(simulation_dir,'t')
    time = force_data[:,4]

    plt.figure(1)
    Ep = get_parameter(simulation_dir,'Ep')
    nup = get_parameter(simulation_dir,'nup')
    E0 = Ep/(1-nup**2)/2
    R0 = 0.03/2
    h = force_data[:, 0]
    F = force_data[:,1]
    plt.plot(h, F)
    plt.ylabel("Force [N]")
    plt.xlabel("Overlap [m]")
    plt.tight_layout()

    plt.figure(2)
    plt.plot(time,F)
    plt.ylabel("Force [N]")
    plt.xlabel("Time [s]")
    plt.tight_layout()

    plt.figure(3)
    plt.plot(time,h)
    plt.ylabel("Overlap [m]")
    plt.xlabel("Time [s]")
    plt.tight_layout()
    plt.show()


if __name__ == '__main__':
    main()