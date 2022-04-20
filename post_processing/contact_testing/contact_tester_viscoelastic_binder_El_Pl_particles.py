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
    force_data = np.genfromtxt('C:/Users/Axel/Documents/DEM/results/contact_testing/plasticity_in_particle/contact_testing/Replicate_EO/contact_testing.dou',
                                delimiter=',')
    simulation_dir = 'C:/Users/Axel/Documents/DEM/DEMsim/simulations/force_model_impact_on_electrode/plastic_particle_contact_test.sim'


    time_parameter = get_parameter(simulation_dir,'t')
    R_eff = get_parameter(simulation_dir,'R')/2
    Syp = get_parameter(simulation_dir,'particle_yield_stress_')
    time = force_data[:,4]
    ticks = time/time_parameter
    plt.figure(1)
    Ep = get_parameter(simulation_dir,'Ep')
    nup = get_parameter(simulation_dir,'nup')
    E0 = Ep/(1-nup**2)/2
    R0 = 0.03/2
    h = force_data[:, 0]
    F = force_data[:,1]
    h_norm = h/(R_eff)
    F_norm = F/(R_eff**2 * Syp)
    plt.plot(h, F)
    plt.ylabel("Force [N]")
    plt.xlabel("Overlap [m]")
    plt.tight_layout()

    # print(ticks)
    plt.figure(2)
    plt.plot(ticks,F)
    plt.ylabel("Force [N]")
    plt.xlabel("Ticks [-]")
    plt.tight_layout()


    plt.figure(3)
    plt.plot(time,F)
    plt.ylabel("Force [N]")
    plt.xlabel("Time [s]")
    plt.tight_layout()

    plt.figure(4)
    plt.plot(time,h)
    plt.ylabel("Overlap [m]")
    plt.xlabel("Time [s]")
    plt.tight_layout()


    plt.figure(5)
    plt.plot(h_norm, F_norm)
    plt.ylabel("Force/R²_eff [-]")
    plt.xlabel("Overlap/R_eff  [-]")
    plt.tight_layout()

    plt.show()

if __name__ == '__main__':
    main()