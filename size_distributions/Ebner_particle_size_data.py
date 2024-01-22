import numpy as np
import pylab as pl
import scipy.special as sps
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import matplotlib
from scipy import integrate, stats, optimize
import pandas as pd


if __name__ == '__main__':
    plt.style.use('axel_style')

    NMC622_Experimental_data_df = pd.read_csv("G:/My Drive/Skola/KTH/PhD/Experiments/Particle size and shape measurement/Q3/q3_lg_EQPC.csv")
    NMC622_q3 = NMC622_Experimental_data_df['q3 lg'].to_numpy()
    NMC622_q3_particle_size = NMC622_Experimental_data_df['xm / μm'].to_numpy()

    NMC622_Experimental_data_df = pd.read_csv("G:/My Drive/Skola/KTH/PhD/Experiments/Particle size and shape measurement/Q3/Q3_EQPC.csv")
    NMC622_Q3 = NMC622_Experimental_data_df['Q3 / %'].to_numpy()
    NMC622_Q3_particle_size = NMC622_Experimental_data_df['xo / μm'].to_numpy()

    NMC622_Experimental_data_df = pd.read_csv(
        "G:/My Drive/Skola/KTH/PhD/Experiments/Particle size and shape measurement/Q0/q0_lg_EQPC.csv")
    NMC622_q0_particle_size = NMC622_Experimental_data_df['xm / μm'].to_numpy()
    NMC622_q0= NMC622_Experimental_data_df['q0 lg'].to_numpy()

    NMC622_Experimental_data_df = pd.read_csv(
        "G:/My Drive/Skola/KTH/PhD/Experiments/Particle size and shape measurement/Q0/Q0_EQPC.csv")
    NMC622_Q0_particle_size = NMC622_Experimental_data_df['xo / μm'].to_numpy()
    NMC622_Q0 = NMC622_Experimental_data_df['Q0 / %'].to_numpy()

    Ebner_particle_raw_data_96wt_df = pd.read_csv(
    "G:/My Drive/Skola/KTH/PhD/Litteratur/X-Ray Tomography of Porous, Transition Metal Oxide Based Lithium Ion Battery Electrodes/NMC_96wt_0bar.dat")
    Ebner_particle_volume_96wt=Ebner_particle_raw_data_96wt_df['volume'].to_numpy() * .37**3
    Ebner_particle_size_96wt = 2*(Ebner_particle_volume_96wt*3/(4*3.14))**(1/3)
    Ebner_particle_size_96wt = Ebner_particle_size_96wt[Ebner_particle_size_96wt<60] #remove all particles larger then 60µm

    Ebner_particle_raw_data_90wt_df = pd.read_csv(
    "G:/My Drive/Skola/KTH/PhD/Litteratur/X-Ray Tomography of Porous, Transition Metal Oxide Based Lithium Ion Battery Electrodes/NMC_90wt_0bar.dat")
    Ebner_particle_volume_90wt=Ebner_particle_raw_data_90wt_df['volume'].to_numpy() * .37**3
    Ebner_particle_size_90wt = 2*(Ebner_particle_volume_90wt*3/(4*3.14))**(1/3)
    Ebner_particle_size_90wt = Ebner_particle_size_90wt[Ebner_particle_size_90wt<60] #remove all particles larger then 60µm

    Laser_diffraction_Experimental_data_from_Ebner_df = pd.read_csv(
        "G:/My Drive/Skola/KTH/PhD/Litteratur/X-Ray Tomography of Porous, Transition Metal Oxide Based Lithium Ion Battery Electrodes/Laser_diffraction_Desnsity_distibution_q3_um-1_Particle_diameter_um.txt")

    particle_diameter_exponent_experimental_laser_diffraction = np.log10(
        Laser_diffraction_Experimental_data_from_Ebner_df['Particle diameter [um]'].to_numpy())
    q_3_experimental_laser_diffraction = Laser_diffraction_Experimental_data_from_Ebner_df[' Density distribution q3 [um^-1]'].to_numpy()

    Tomography_Experimental_data_from_Ebner_df = pd.read_csv(
        "G:/My Drive/Skola/KTH/PhD/Litteratur/X-Ray Tomography of Porous, Transition Metal Oxide Based Lithium Ion Battery Electrodes/Tomography_Desnsity_distibution_q3_um-1_Particle_diameter_um.txt")

    particle_diameter_exponent_experimental_tomography = np.log10(
        Tomography_Experimental_data_from_Ebner_df['Particle diameter [um]'].to_numpy())
    q_3_experimental_tomography = Tomography_Experimental_data_from_Ebner_df[
        ' Density distribution q3 [um^-1]'].to_numpy()

    # =NMC622 log CDF/PDF for Q0 an Q3 =====================================================================================
    fig_q0_q3, ax_q0_q3  = plt.subplots()
    lns_q3_NMC622 = ax_q0_q3.plot((NMC622_q3_particle_size), NMC622_q3, label='NMC622 - q3')
    lns_q0_NMC622 = ax_q0_q3.plot((NMC622_q0_particle_size), NMC622_q0, label='NMC622 - q0')
    ax_q0_q3.set_xscale('log')

    ax_q0_q3.xaxis.set_major_formatter(matplotlib.ticker.LogFormatter())
    ax_q0_q3.xaxis.set_minor_formatter(matplotlib.ticker.ScalarFormatter())
    ax_q0_q3.set_xlim(xmin=1,xmax=60)
    ax_q0_q3.set_ylim(ymin=0)
    ax_q0_q3.set_xlabel('log (particle diameter) [µm]')
    ax_q0_q3.set_ylabel('Density distribution $q_i$ $[µm^-1]$')
    # ax_q0_q3.legend(loc='best')
    ax_Q0_Q3 = ax_q0_q3.twinx()
    lns_Q3_NMC622 = ax_Q0_Q3.plot((NMC622_Q3_particle_size), NMC622_Q3,'--', label='NMC622 - Q3')
    lns_Q0_NMC622 = ax_Q0_Q3.plot((NMC622_Q0_particle_size), NMC622_Q0,'--', label='NMC622 - Q0')
    ax_Q0_Q3.set_ylabel('Cumulative distribution $Q_i$ $[µm^{-1}]$')
    ax_Q0_Q3.set_ylim(ymin=0,ymax=100)
    lns_Qq = lns_Q3_NMC622 + lns_q3_NMC622+ lns_Q0_NMC622 + lns_q0_NMC622  # + line_q_0_fitted
    labels_Qq = [l.get_label() for l in lns_Qq]
    ax_q0_q3.legend(lns_Qq, labels_Qq, loc='upper left')

    # =LOGARITHMIC PDF q3 FROM EBNER==================================================================================

    fig_q3_xlog, ax_q3_xlog = plt.subplots()
    lns_q3_laser_diffraction_log = ax_q3_xlog.plot(10** particle_diameter_exponent_experimental_laser_diffraction,
                                                  q_3_experimental_laser_diffraction, label='Laser diffraction - q3')
    lns_q3_tomography_log = ax_q3_xlog.plot(10**particle_diameter_exponent_experimental_tomography,
                                                  q_3_experimental_tomography, label='Tomography - q3')
    lns_q3_NMC622_log = ax_q3_xlog.plot((NMC622_q3_particle_size), NMC622_q3, label = 'NMC622 - q3')
    ax_q3_xlog.set_xscale('log')
    ax_q3_xlog.set_xlabel('log (particle diameter) [µm]')
    ax_q3_xlog.set_ylabel('Density distribution $q_3$ $[µm^-1]$')
    ax_q3_xlog.legend(loc='best')

    # =LOG PDF q0=======================================================================================================
    fig_q0_xlog, ax_q0_xlog = plt.subplots()
    # lns_q3_laser_diffraction_log = ax_q3_xlog.plot(10** particle_diameter_exponent_experimental_laser_diffraction,
    #                                               q_3_experimental_laser_diffraction, label='Laser diffraction - q3')
    # lns_q3_tomography_log = ax_q3_xlog.plot(10**particle_diameter_exponent_experimental_tomography,
    #                                               q_3_experimental_tomography, label='Tomography - q3')
    lns_q0_NMC622_log = ax_q0_xlog.plot((NMC622_q0_particle_size), NMC622_q0, label = 'NMC622 - q0')
    lns_q3_NMC622_log = ax_q0_xlog.plot((NMC622_q3_particle_size), NMC622_q3, label = 'NMC622 - q3')

    ax_q0_xlog.set_xscale('log')
    ax_q0_xlog.set_xlabel('log (particle diameter) [µm]')
    ax_q0_xlog.set_ylabel('Density distribution $q_0$ $[µm^-1]$')
    ax_q0_xlog.legend(loc='best')

    # # =LOGARITHMIC PDF q3 FROM EBNER==================================================================================
    # fig_q3_log, ax_q3_log = plt.subplots()
    # lns_q3_laser_diffraction_log = ax_q3_log.plot(particle_diameter_exponent_experimental_laser_diffraction,
    #                                               q_3_experimental_laser_diffraction, label='Laser diffraction - q3')
    # lns_q3_tomography_log = ax_q3_log.plot(particle_diameter_exponent_experimental_tomography,
    #                                               q_3_experimental_tomography, label='Tomography - q3')
    # lns_q3_NMC622_log = ax_q3_log.plot(np.log10(NMC622_q3_particle_size), NMC622_q3, label = 'NMC622 - q3')
    # ax_q3_log.set_xlabel('log (particle diameter) [µm]')
    # ax_q3_log.set_ylabel('Density distribution $q_3$ $[µm^-1]$')
    # ax_q3_log.legend(loc='best')


    # =rescaleing by http: // www.ebyte.it / library / docs / math04a / PdfChangeOfCoordinates04.html===================
    # pdf_lin = d_pdf. / (log(10) * 10. ^ d_exponent);
    q3_lin_experimental_laser_diffraction = q_3_experimental_laser_diffraction / \
                                            (np.log(10) * 10 ** particle_diameter_exponent_experimental_laser_diffraction)
    q_3_lin_experimental_tomography = q_3_experimental_tomography / \
                                      (np.log(10) * 10 ** particle_diameter_exponent_experimental_tomography)
    q_3_lin_experimental_NMC622 = NMC622_q3 / \
                                      (np.log(10) * NMC622_q3_particle_size)
    q_0_lin_experimental_NMC622 = NMC622_q0 / \
                                  (np.log(10) * NMC622_q0_particle_size)

    print(integrate.cumtrapz(q_3_lin_experimental_NMC622, (NMC622_q3_particle_size)))




    (Ebner_particle_size_96wt_shape,Ebner_particle_size_96wt_loc,Ebner_particle_size_96wt_scale) = stats.gamma.fit(Ebner_particle_size_96wt,loc=min(Ebner_particle_size_96wt))
    Ebner_size_array = np.linspace(1,60,10000)
    Ebner_96wt_q0_pdf = stats.gamma.pdf(Ebner_size_array, Ebner_particle_size_96wt_shape, Ebner_particle_size_96wt_loc, Ebner_particle_size_96wt_scale)
    Ebner_96wt_q3_pdf = Ebner_size_array ** 3 * Ebner_96wt_q0_pdf * len(Ebner_particle_size_96wt) / np.sum(Ebner_particle_size_96wt ** 3)


    # =NMC622 lin CDF/PDF for Q0 an Q3 =====================================================================================
    fig_lin_q0_q3, ax_lin_q0_q3  = plt.subplots()
    lns_lin_q3_NMC622 = ax_lin_q0_q3.plot((NMC622_q3_particle_size), q_3_lin_experimental_NMC622, label='NMC622 - q3')
    lns_lin_q0_NMC622 = ax_lin_q0_q3.plot((NMC622_q0_particle_size), q_0_lin_experimental_NMC622, label='NMC622 - q0')
    # ax_q0_q3.set_xscale('log')

    # ax_q0_q3.xaxis.set_major_formatter(matplotlib.ticker.LogFormatter())
    # ax_q0_q3.xaxis.set_minor_formatter(matplotlib.ticker.ScalarFormatter())
    ax_lin_q0_q3.set_xlim(xmin=1,xmax=30)
    ax_lin_q0_q3.set_ylim(ymin=0)
    ax_lin_q0_q3.set_xlabel('particle diameter [µm]')
    ax_lin_q0_q3.set_ylabel('Density distribution $q_i$ $[µm^-1]$')
    # ax_q0_q3.legend(loc='best')
    ax_lin_Q0_Q3 = ax_lin_q0_q3.twinx()
    lns_lin_Q3_NMC622 = ax_lin_Q0_Q3.plot((NMC622_Q3_particle_size), NMC622_Q3,'--', label='NMC622 - Q3')
    lns_lin_Q0_NMC622 = ax_lin_Q0_Q3.plot((NMC622_Q0_particle_size), NMC622_Q0,'--', label='NMC622 - Q0')
    ax_lin_Q0_Q3.set_ylabel('Cumulative distribution $Q_i$ $[µm^{-1}]$')
    ax_lin_Q0_Q3.set_ylim(ymin=0,ymax=100)
    lns_Qq = lns_lin_Q3_NMC622 + lns_lin_q3_NMC622+ lns_lin_Q0_NMC622 + lns_lin_q0_NMC622  # + line_q_0_fitted
    labels_Qq = [l.get_label() for l in lns_Qq]
    ax_lin_q0_q3.legend(lns_Qq, labels_Qq, loc='upper left')

    # =FIG HISTOGRAM EBNER PARTICLE DATA================================================================================
    fig_particle_histogram,ax_particle_histogram = plt.subplots()
    lns_Ebner = ax_particle_histogram.hist([Ebner_particle_size_96wt,Ebner_particle_size_90wt], 100, density=True, label=['$Tomography, raw data - 96wt$','$Tomography, raw data - 90wt$'])
    lns_q0_ebner = ax_particle_histogram.plot(Ebner_size_array,
                                   Ebner_96wt_q0_pdf, label='Tomography, raw data- 96wt - fitted q0')
    # lns_Ebner2 = ax_particle_histogram.hist(Ebner_particle_size2, 100,align='right',rwidth=.5, density=True, label='$Ebner2$')

    lnd_q0_NMC622 = ax_particle_histogram.plot(NMC622_q0_particle_size,q_0_lin_experimental_NMC622, label='NMC622 q0')

    ax_particle_histogram.set_title('Particle size distribution from tomography raw data')
    ax_particle_histogram.set_xlabel('particle diameter [µm]')
    ax_particle_histogram.set_ylabel('Number of particles normalised/ q0')
    ax_particle_histogram.legend(loc='best')


    # =FIG q3 CDF FROM FIGURE AND FITTED DATA===========================================================================
    fig_q3, ax_q3 = plt.subplots()
    lns_q3_laser_diffraction = ax_q3.plot(10**particle_diameter_exponent_experimental_laser_diffraction,
                                                  q3_lin_experimental_laser_diffraction, label='Laser diffraction - q3')
    lns_q3_tomography = ax_q3.plot(10**particle_diameter_exponent_experimental_tomography,
                                           q_3_lin_experimental_tomography, label='Tomography - q3')
    lns_q3_NMC622 = ax_q3.plot(NMC622_q3_particle_size,
                                           q_3_lin_experimental_NMC622, label='NMC622 - q3')



    lns_q3_ebner = ax_q3.plot(Ebner_size_array,
                                   Ebner_96wt_q3_pdf, label='Tomography, raw data - fitted q3')
    ax_q3.set_xlabel('particle diameter [µm]')
    ax_q3.set_ylabel('Density distribution $q_i$ $[µm^{-1}]$')
    ax_q3.legend(loc='best')
    plt.show()