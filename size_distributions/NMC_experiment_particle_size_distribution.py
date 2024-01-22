# ==========================
# Completed by Axel 23/10/06
# ==========================

import numpy as np
import pylab as pl
import scipy.special as sps
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import matplotlib
from scipy import integrate, stats, optimize
import pandas as pd
import random


class Point:
    def __init__(self, x, y):
        self.x = x
        self.y = y


def pdf_q0(xp,fp,x):
    return np.interp(x, xp, fp)

def rand_point():
    rand_x = random.uniform(0, 1)
    rand_y = random.uniform(0, 1)
    return (rand_x,rand_y)

if __name__ == '__main__':
    plt.style.use('axel_style')

    NMC622_q0_Experimental_data_df = pd.read_csv(
        "G:/My Drive/Skola/KTH/PhD/Experiments/Particle size and shape measurement/Q0/q0_lg_EQPC.csv")
    NMC622_q0_particle_size = NMC622_q0_Experimental_data_df['xm / μm'].to_numpy()
    NMC622_q0 = NMC622_q0_Experimental_data_df['q0 lg'].to_numpy()
    NMC622_q0_lin = NMC622_q0 / (np.log(10) * NMC622_q0_particle_size)

    NMC622_Q0_Experimental_data_df = pd.read_csv(
        "G:/My Drive/Skola/KTH/PhD/Experiments/Particle size and shape measurement/Q0/Q0_EQPC.csv")
    NMC622_Q0_particle_size = NMC622_Q0_Experimental_data_df['xo / μm'].to_numpy()
    NMC622_Q0 = NMC622_Q0_Experimental_data_df['Q0 / %'].to_numpy()

    NMC622_Experimental_data_df = pd.read_csv("G:/My Drive/Skola/KTH/PhD/Experiments/"
                                              "Particle size and shape measurement/Q3/q3_lg_EQPC.csv")
    NMC622_q3 = NMC622_Experimental_data_df['q3 lg'].to_numpy()
    NMC622_q3_particle_size = NMC622_Experimental_data_df['xm / μm'].to_numpy()
    NMC622_q3_lin = NMC622_q3 / (np.log(10) * NMC622_q3_particle_size)

    NMC622_Experimental_data_df = pd.read_csv("G:/My Drive/Skola/KTH/PhD/Experiments/"
                                              "Particle size and shape measurement/Q3/Q3_EQPC.csv")
    NMC622_Q3 = NMC622_Experimental_data_df['Q3 / %'].to_numpy()
    NMC622_Q3_particle_size = NMC622_Experimental_data_df['xo / μm'].to_numpy()

    # =PLOT OF q0, Q0 CALCULATED AND Q0 FROM EXPERIMENTS============================================================
    fig_q0_Q0, ax_q0_Q0  = plt.subplots()
    lns_q0 = ax_q0_Q0.plot(NMC622_q0_particle_size, NMC622_q0,'-*', label='q0')
    ax_q0_Q0.set_xscale('log')

    ax_q0_Q0.xaxis.set_major_formatter(matplotlib.ticker.LogFormatter())
    ax_q0_Q0.xaxis.set_minor_formatter(matplotlib.ticker.ScalarFormatter())
    ax_q0_Q0.set_xlim(xmin=1,xmax=60)
    ax_q0_Q0.set_ylim(ymin=0)
    ax_q0_Q0.set_xlabel('log (particle diameter) [µm]')
    ax_q0_Q0.set_ylabel('Density distribution $q_i$ $[µm^-1]$')
    # ax_q0_q3.legend(loc='best')
    ax_Q0 = ax_q0_Q0.twinx()
    ax_Q0.set_xscale('log')
    lns_Q0_calc = ax_Q0.plot(NMC622_q0_particle_size[1:],
                             integrate.cumtrapz(NMC622_q0, np.log10(NMC622_q0_particle_size))*100,
                             '--', label='Q0 - calc')
    lns_Q0 = ax_Q0.plot((NMC622_Q0_particle_size), NMC622_Q0,'--', label='Q0')
    ax_Q0.set_ylabel('Cumulative distribution $Q_i$ $[µm^{-1}]$')
    ax_Q0.set_ylim(ymin=0,ymax=100)
    lns_Qq = lns_q0 + lns_Q0_calc + lns_Q0
    labels_Qq = [l.get_label() for l in lns_Qq]
    ax_q0_Q0.legend(lns_Qq, labels_Qq, loc='center right')
    ax_Q0.xaxis.set_major_formatter(matplotlib.ticker.LogFormatter())
    ax_Q0.xaxis.set_minor_formatter(matplotlib.ticker.ScalarFormatter())
    # ==================================================================================================================

    # Create points in xy plane
    n_points = 10**5
    x_point  = []
    y_point = []
    for i in range(n_points):
       x_point.append(min(NMC622_q0_particle_size) + (max(NMC622_q0_particle_size) - min(NMC622_q0_particle_size)) *
                      rand_point()[0])
       y_point.append(max(NMC622_q0_lin) * rand_point()[1])
    x_point = np.asarray(x_point)
    y_point = np.asarray(y_point)

    x = np.linspace(min(NMC622_q0_particle_size),max(NMC622_q0_particle_size),100)
    # plt.plot(x,pdf_q0(NMC622_q0_particle_size,NMC622_q0_lin,x),'-*')
    # plt.plot(NMC622_q0_particle_size,NMC622_q0_lin,'-x')

    # sort out points under q0
    l = []
    for i in range(len(x_point)):
        if y_point[i] <= pdf_q0(NMC622_q0_particle_size,NMC622_q0_lin,x_point[i]):
            l.append(x_point[i])
    valid_points = np.asarray(l)
    counts, bins = np.histogram(valid_points,bins=100,density=True)

    # =PDF AND HISTOGRAM FOR q0=========================================================================================
    fig_q0, ax_q0 = plt.subplots()
    # lns_pdf = ax_q0.plot(x,pdf_q0(NMC622_q0_particle_size,NMC622_q0_lin,x))
    lns_pdf = ax_q0.plot(NMC622_q0_particle_size,NMC622_q0_lin,'-*')
    # ax_q0.scatter(x_point, y_point,alpha=0.1)
    ax_q0.stairs(counts, bins,fill=True)
    ax_q0.set_xlabel('Particle size [µm]')
    ax_q0.set_ylabel('Probability density $q_0$ $[µm^{-1}]$')

    # Calculations for volume histogram
    counts_vol = np.zeros(len(counts))
    for i in valid_points:
        j = 0
        if i >= bins[-1]:
            j = len(bins)-2
        else:
            while i >= bins[j]:
                j += 1
            j -= 1
        counts_vol[j] += i**3
    counts_vol = (1/(bins[1]-bins[0])) * counts_vol / np.sum(valid_points**3)

    # =PDF AND HISTOGRAM FOR q3=========================================================================================
    fig_q3, ax_q3 = plt.subplots()
    lns_pdf = ax_q3.plot(NMC622_q3_particle_size,NMC622_q3_lin,'-*')
    # lns_pdf = ax_q3.plot(x, pdf_q0(NMC622_q3_particle_size, NMC622_q3_lin, x))
    ax_q3.stairs(counts_vol, bins,fill=True)
    # print(integrate.cumtrapz(NMC622_q3_lin, (NMC622_q3_particle_size)))
    ax_q3.set_xlabel('Particle size [µm]')
    ax_q3.set_ylabel('Probability density $q_3$ $[µm^{-1}]$')

    # =CREATION OF PARTICLE SIZE FILE===================================================================================
    d_min = 9.10
    d_max = 23.2
    scale_factor = 1/100
    diameter_to_radius = 1/2
    write_points = valid_points[valid_points < d_max]
    write_points = write_points[write_points > d_min] * scale_factor * diameter_to_radius
    ###np.savetxt("G:/My Drive/Skola/KTH/PhD/Experiments/Particle size and shape measurement/Particle_size.dat",write_points,fmt='%.8f')
    # =END CREATION OF PARTICLE SIZE FILE===============================================================================

    # =COMPARE SIZE DISTRIBUTIONS=======================================================================================
    old_particle_file = 'C:/Users/Axel/Documents/DEM/DEMsim/simulations/article_1_final_runs/particle_radius_distr_gamma_k_492_theta_263_r_min_4_r_max_12_cm_size.dat'
    old_particle_data = np.genfromtxt(old_particle_file)

    new_particle_file = "G:/My Drive/Skola/KTH/PhD/Experiments/Particle size and shape measurement/Particle_size.dat"
    new_particle_data = np.genfromtxt(new_particle_file)

    fig_old_new, ax_old_new = plt.subplots()
    ax_old_new.hist(new_particle_data,bins=100,density=True,alpha=0.99,label='new size distribution')
    ax_old_new.hist(old_particle_data,bins=100,density=True,alpha=0.5,label='old size distribution')
    ax_old_new.legend(loc='best')
    ax_old_new.set_xlabel('Particle radius [m]')
    ax_old_new.set_ylabel('Density')

    print('show plots')
    plt.show()