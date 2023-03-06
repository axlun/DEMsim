import numpy as np
import pylab as pl
import scipy.special as sps
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import matplotlib
from scipy import integrate, stats, optimize
import pandas as pd


# plt.rcParams['text.usetex'] = True

# Particles size disribution of NMC particles
# Gamma distribution with alpha=4.4152, beta=1.9944, gives particle diameter in m

# Gamma distribution with Shape  k=4.9242 and scale theta=2.6274, gives particle diameter in m

def f_gamma_pdf(x, kappa, theta):
    return x ** (kappa - 1) * (np.exp(-x / theta) / (sps.gamma(kappa) * theta ** kappa))


def f_gamma_cdf(x, kappa, theta):
    return sps.gammainc(kappa, (x / theta))  # (1/Gamma(kappa)) * gamma(kappa,x/theta)


if __name__ == '__main__':
    plt.style.use('axel_style')

    Ebner_particle_raw_data_df = pd.read_csv(
    "G:/My Drive/Skola/KTH/PhD/Litteratur/X-Ray Tomography of Porous, Transition Metal Oxide Based Lithium Ion Battery Electrodes/NMC_96wt_0bar.dat")
    Ebner_particle_volume=Ebner_particle_raw_data_df['volume'].to_numpy()
    Ebner_particle_size = 2*(Ebner_particle_volume*3/(4*3.14))**(1/3)
    Ebner_particle_size = Ebner_particle_size[Ebner_particle_size<60] #remove all particles larger then 60µm

    Ebner_particle_raw_data_df2 = pd.read_csv(
    "G:/My Drive/Skola/KTH/PhD/Litteratur/X-Ray Tomography of Porous, Transition Metal Oxide Based Lithium Ion Battery Electrodes/NMC_90wt_0bar.dat")
    Ebner_particle_volume2=Ebner_particle_raw_data_df2['volume'].to_numpy()
    Ebner_particle_size2 = 2*(Ebner_particle_volume2*3/(4*3.14))**(1/3)
    Ebner_particle_size2 = Ebner_particle_size2[Ebner_particle_size2<60] #remove all particles larger then 60µm

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

    # =LOGARITHMIC PDF q3 FROM EBNER==================================================================================
    fig_q3_log, ax_q3_log = plt.subplots()
    lns_q3_laser_diffraction_log = ax_q3_log.plot(particle_diameter_exponent_experimental_laser_diffraction,
                                                  q_3_experimental_laser_diffraction, label='Laser diffraction - q3')
    lns_q3_tomography_log = ax_q3_log.plot(particle_diameter_exponent_experimental_tomography,
                                                  q_3_experimental_tomography, label='Tomography - q3')
    ax_q3_log.set_xlabel('log (particle diameter) [µm]')
    ax_q3_log.set_ylabel('Density distribution $q_3$ $[µm^-1]$')

    # =rescaleing by http: // www.ebyte.it / library / docs / math04a / PdfChangeOfCoordinates04.html===================
    # pdf_lin = d_pdf. / (log(10) * 10. ^ d_exponent);
    q3_lin_experimental_laser_diffraction = q_3_experimental_laser_diffraction / \
                                            (np.log(10) * 10 ** particle_diameter_exponent_experimental_laser_diffraction)
    q_3_lin_experimental_tomography = q_3_experimental_tomography / \
                                      (np.log(10) * 10 ** particle_diameter_exponent_experimental_tomography)

    # =FITTING EBNER DATA TO PDF/CDF====================================================================================
    (Ebner_shape,Ebner_loc,Ebner_scale) = stats.gamma.fit(Ebner_particle_size,loc=min(Ebner_particle_size))
    print("Ebner fit (shape, loc, scale): " + str((Ebner_shape,Ebner_loc,Ebner_scale)))
    print('Ebner min particle = ' + str(min(Ebner_particle_size)))
    Ebner_size_array = np.linspace(1,60,10000)
    Ebner_q0_pdf = stats.gamma.pdf(Ebner_size_array, Ebner_shape, Ebner_loc, Ebner_scale)
    Ebner_q3_pdf = Ebner_size_array ** 3 * Ebner_q0_pdf * len(Ebner_particle_size) / np.sum(Ebner_particle_size ** 3)
    Ebner_q_0_sample = stats.gamma.rvs(Ebner_shape, Ebner_loc, Ebner_scale, size=10000)
    # print(Ebner_size_array)
    Ebner_q0_pdf_log = Ebner_size_array * Ebner_q0_pdf/np.log10(2.718)
    Ebner_size_array_log = np.log10(Ebner_size_array)
    lns_q0_Ebner_log = ax_q3_log.plot(Ebner_size_array_log,
                                                  Ebner_q0_pdf_log, label='Tomography, raw data - fitted q0')
    ax_q3_log.legend(loc='best')

    # =FIG HISTOGRAM EBNER PARTICLE DATA================================================================================
    fig_particle_histogram,ax_particle_histogram = plt.subplots()
    lns_Ebner = ax_particle_histogram.hist([Ebner_particle_size,Ebner_particle_size2], 100, density=True, label=['$Tomography, raw data - 96wt$','$Tomography, raw data - 90wt$'])
    lns_q0_ebner = ax_particle_histogram.plot(Ebner_size_array,
                                   Ebner_q0_pdf, label='Tomography, raw data - fitted q0')
    # lns_Ebner2 = ax_particle_histogram.hist(Ebner_particle_size2, 100,align='right',rwidth=.5, density=True, label='$Ebner2$')
    ax_particle_histogram.set_title('Particle size distribution from tomography raw data')
    ax_particle_histogram.set_xlabel('particle diameter [µm]')
    ax_particle_histogram.set_ylabel('Number of particles normalised/ q0')
    ax_particle_histogram.legend(loc='best')

    # =FIG q3 LINEAR===================================================================================================
    fig_q3, ax_q3 = plt.subplots()
    lns_q3_laser_diffraction = ax_q3.plot(10**particle_diameter_exponent_experimental_laser_diffraction,
                                                  q3_lin_experimental_laser_diffraction, label='Laser diffraction - q3')
    lns_q3_tomography = ax_q3.plot(10**particle_diameter_exponent_experimental_tomography,
                                           q_3_lin_experimental_tomography, label='Tomography - q3')
    lns_q0_ebner = ax_q3.plot(Ebner_size_array,
                                   Ebner_q0_pdf, label='Tomography, raw data - fitted q0')
    ax_q3.set_xlabel('particle diameter [µm]')
    ax_q3.set_ylabel('Density distribution $q_i$ $[µm^{-1}]$')
    ax_q3.legend(loc='best')
    print('Integration of q3_lin_experimental_laser_diffraction = ' + str(integrate.cumtrapz(q3_lin_experimental_laser_diffraction,10**particle_diameter_exponent_experimental_laser_diffraction)[-1]))
    print('Integration of q3_lin_experimental_tomography = ' + str(integrate.cumtrapz(q_3_lin_experimental_tomography,10**particle_diameter_exponent_experimental_tomography)[-1]))

    # shape_q3_guess = 5
    # loc_q3_guess = 0
    # scale_q3_guess = 3
    #
    # q_3_guess = stats.gamma.pdf(10 ** particle_diameter_exponent_experimental_laser_diffraction, shape_q3_guess, loc_q3_guess, scale_q3_guess)
    #
    # # q_3_guess = f_gamma_pdf(10 ** particle_diameter_exponent_experimental_laser_diffraction, shape_q3_guess, scale_q3_guess)
    #
    # [shape_q_3_opt,  loc_q_3_opt,scale_q_3_opt], pcov = optimize.curve_fit(
    #     f=stats.gamma.pdf,  # model function
    #     xdata=10 ** particle_diameter_exponent_experimental_laser_diffraction,  # x data
    #     ydata=q3_lin_experimental_laser_diffraction,  # y data
    #     p0=(shape_q3_guess, loc_q3_guess, scale_q3_guess),  # initial value of the parameters
    # )
    # # q_3_fitted = f_gamma_pdf(10 ** particle_diameter_exponent_experimental_laser_diffraction, shape_q_3_opt, scale_q_3_opt)
    # q_3_fitted = stats.gamma.pdf(10 ** particle_diameter_exponent_experimental_laser_diffraction, shape_q_3_opt,
    #                              loc_q_3_opt, scale_q_3_opt)
    #
    # # =EXPERIMENTAL PDF q3==============================================================================================
    # fig_q_3, ax_q_3 = plt.subplots()
    # q_3_experimental_laser_diffraction = ax_q_3.plot(10 ** particle_diameter_exponent_experimental_laser_diffraction, q3_lin_experimental_laser_diffraction, label='Experimental values')
    # q_3_guess = ax_q_3.plot(10 ** particle_diameter_exponent_experimental_laser_diffraction, q_3_guess, label='Guess')
    # q_3_fit = ax_q_3.plot(10 ** particle_diameter_exponent_experimental_laser_diffraction, q_3_fitted, label='Fitted')
    # plt.legend(loc='best')
    # print('Gamma parameters for q3 fitted to experiment\n' + 'kappa_q_3 = ' + str(shape_q_3_opt) + '\nloc_q_3 = '
    #       + str(loc_q_3_opt) + '\nscale_q_3 = ' + str(scale_q_3_opt) + '\n')

    # D_min = 1.5
    # D_max = 60
    # D_array = np.linspace(D_min, D_max, num=5000)
    # q_3 = f_gamma_pdf(D_array, shape_q_3_opt, scale_q_3_opt)
    # Q_3 = np.append(0, integrate.cumtrapz(q_3, D_array))
    #
    # # Eq. 25 in Stoyan and Unland 2022
    # Q_0 = np.array([])
    # for x in D_array[1:]:
    #     Q_0_integration = \
    #     integrate.cumtrapz(D_array[D_array <= x] ** -3 * D_array[D_array <= x] ** (shape_q_3_opt - 1) * (
    #             np.exp(-D_array[D_array <= x] / scale_q_3_opt) / (
    #                 sps.gamma(shape_q_3_opt) * scale_q_3_opt ** shape_q_3_opt)), D_array[D_array <= x])[-1] / \
    #     integrate.cumtrapz(D_array ** -3 * D_array ** (shape_q_3_opt - 1) *
    #                        (np.exp(-D_array / scale_q_3_opt) / (
    #                                    sps.gamma(shape_q_3_opt) * scale_q_3_opt ** shape_q_3_opt)), D_array)[-1]
    #     Q_0 = np.append(Q_0, Q_0_integration)
    # Q_0 = np.append(0, Q_0)
    #
    # q_0 = np.gradient(Q_0, D_array)
    # # print(integrate.cumtrapz(q_0, D_array))
    # shape_q_0_guess = 2.5
    # scale_q_0_guess = 1.5
    #
    # q_0_guess = f_gamma_pdf(D_array, scale_q_0_guess, scale_q_0_guess)
    #
    # [shape_q_0_opt, scale_q_0_opt], pcov = optimize.curve_fit(
    #     f=f_gamma_pdf,  # model function
    #     xdata=D_array,  # x data
    #     ydata=q_0,  # y data
    #     p0=(shape_q_0_guess, scale_q_0_guess),  # initial value of the parameters
    # )
    # # print(popt)
    # print('Gamma parameters for q_0\n' + 'shape_q_0 = ' + str(shape_q_0_opt) + '\nscale_q_0 = '
    #       + str(scale_q_0_opt))
    #
    # # =FIGURE OF CDF AND PDF============================================================================================
    # fig_q, ax_q = plt.subplots()
    # line_q_3 = ax_q.plot(D_array, q_3, label=r'$q_3(D)$')  # ,marker ='*', markevery=50, markersize=8,
    # line_q_0 = ax_q.plot(D_array, q_0, label=r'$q_0(D)$')
    # # line_q_0_fitted = ax_q.plot(D_array,f_gamma_cdf(D_array,shape_q_0_opt,scale_q_0_opt),'--',label=r'$q_0(D)$ fitted')
    # # line_q0_guess = ax_Q.plot(D_array,q_0_guess,label=r'$q_0$ guess')
    # ax_Q = ax_q.twinx()
    # line_Q3 = ax_Q.plot(D_array, Q_3, marker='o', markevery=50, markersize=8, label=r'$Q_3(D)$')
    # line_Q0 = ax_Q.plot(D_array, Q_0, marker='X', markevery=50, markersize=8, label=r'$Q_0(D)$')
    #
    # ax_q.set_xlim(xmin=0, xmax=50)
    # ax_q.set_ylim(ymin=0)  # ,ymax=1.1)
    # ax_q.set_xlabel('Particle diameter D [µm]')
    # ax_q.set_ylabel('Density distribution $q_i$ $[µm^{-1}]$')
    # ax_Q.set_ylim(ymin=0, ymax=1)
    # ax_Q.set_ylabel('Cumulative distribution $Q_i$ $[µm^{-1}]$')
    #
    # lns_Qq = line_Q3 + line_Q0 + line_q_3 + line_q_0  # + line_q_0_fitted
    # labels_Qq = [l.get_label() for l in lns_Qq]
    # ax_q.legend(lns_Qq, labels_Qq, loc='best')
    #
    # s_q_0 = np.random.gamma(shape_q_0_opt, scale_q_0_opt, 18692)
    # s_q_3 = np.random.gamma(shape_q_3_opt, scale_q_3_opt, 18692)
    #
    # V_tot = np.sum(4 / 3 * 3.14 * (.5 * s_q_0) ** 3)
    # V_part = np.sum(4 / 3 * 3.14 * (.5 * s_q_0[s_q_0 < 12.073]) ** 3)
    # # print(V_part / V_tot)
    #
    # bins = np.arange(0,60.1,0.5)
    # # print(bins)
    # # =HISTOGRAM FOR Q0 AND Q3 AND EBNER DIRSTRIBUTIONS=================================================================
    # fig_histogram, ax_histogram = plt.subplots(1,3)
    # lns_q_3 = ax_histogram[0].hist(s_q_3, bins, density=True, label='$q_3$')
    # lns_q_0 = ax_histogram[1].hist(s_q_0, bins, density=True, label='$q_0$')
    # lns_Ebner = ax_histogram[2].hist(Ebner_particle_size, bins, density=True, label='$Ebner$')
    # ax_histogram[0].set_title('q_3')
    # ax_histogram[1].set_title('q_0')
    # ax_histogram[2].set_title('Ebner')
    # for x in ax_histogram:
    #     x.set_xlabel('Particle size [µm]')
    #     x.set_ylabel('Number of particles')
    #     # x.legend(loc='best')
    # # ax_histogram.set_xlabel('Particle diameter [µm]')
    # # ax_histogram.set_ylabel('Nomalised number of particles')
    #
    #
    # # =FIG qo from Ebner================================================================================================
    # fig_ebner_fit_hist, ax_ebner_fit_hist = plt.subplots()
    # lns_Ebner_hist = ax_ebner_fit_hist.hist(Ebner_particle_size, bins, density=True, label='Histogram')
    # lns_Ebner_pdf = ax_ebner_fit_hist.plot(Ebner_size_array, Ebner_q0_pdf, label='PDF')
    # ax_ebner_fit_hist.set_xlabel('Particle size [µm]')
    # ax_ebner_fit_hist.set_ylabel('PDF q0')
    #
    #
    # # =FIG OF ALL PDFS qo/q3/Ebner======================================================================================
    # fig_q0_q3_ebner, ax_q0_q3_ebner = plt.subplots()
    # # lns_q0 = ax_q0_q3_ebner.plot(D_array, q_0,label='q0')
    # # lns_q3 = ax_q0_q3_ebner.plot(D_array, q_3,label='q3')
    # lns_q3_tomography = ax_q0_q3_ebner.plot(10 **particle_diameter_exponent_experimental_tomography, q_3_lin_experimental_tomography, label='Tomography')
    # lns_q3_laser_diffraction = ax_q0_q3_ebner.plot(10 **particle_diameter_exponent_experimental_laser_diffraction, q3_lin_experimental_laser_diffraction, label='Laser diffraction')
    # lns_q3_Ebner = ax_q0_q3_ebner.plot(Ebner_size_array, Ebner_q3_pdf, label='Ebner q3')
    #
    #
    # lns_Ebner = ax_q0_q3_ebner.plot(Ebner_size_array, Ebner_q0_pdf, label='Ebner')
    # ax_q0_q3_ebner.legend(loc='best')
    # ax_q0_q3_ebner.set_xlabel('Particle size [µm]')
    # ax_q0_q3_ebner.set_ylabel('PDF')
    #
    #
    # # =HISTOGRAM FOR Q0 AND Q3 DIRSTRIBUTIONS===========================================================================
    # fig_histogram, ax_histogram = plt.subplots()
    # lns_q_3 = ax_histogram.hist(s_q_3, 100, density=True, label='$q_3$')
    # lns_q_0 = ax_histogram.hist(s_q_0, 100, density=True, label='$q_0$')
    # ax_histogram.set_xlabel('Particle diameter [µm]')
    # ax_histogram.set_ylabel('Nomalised number of particles')
    # ax_histogram.legend(loc='best')
    #
    # D_min = 6
    # D_max = 20
    # radius_data_for_simulation = 1E-2 * s_q_0[(s_q_0 > D_min) & (s_q_0 < D_max)] / 2
    #
    # V_part_reduced_radius = np.sum(
    #     4 / 3 * 3.14 * ((1E2 * radius_data_for_simulation[radius_data_for_simulation < (12.073 / 2)]) ** 3))

    # print(V_part_reduced_radius / V_tot)
    #
    # print(max(radius_data_for_simulation))
    # print(min(radius_data_for_simulation))
    # np.savetxt("G:/My Drive/Skola/KTH/PhD/Litteratur/X-Ray Tomography of Porous, Transition Metal Oxide Based"
               # " Lithium Ion Battery Electrodes/Particle_radius_data_gamma_k_162_theta_267_r_min_3_r_max_10_from_q_0_distribution.dat",radius_data_for_simulation,fmt='%.8f')

    # fig_histogram_particle_radius, ax_histogram_particle_radius = plt.subplots()
    # lns_particle_radius = ax_histogram_particle_radius.hist(radius_data_for_simulation, 50)
    # ax_histogram_particle_radius.set_title('Histogram for particle radius distribution')
    # ax_histogram_particle_radius.set_xlabel('Particle radius [µm]')

    # Distrubtion from purchased NMC622
    Q_0_NMC622_experimental = np.array([0.01, 10, 50, 90, 100]) / 100
    D_NMC622_experimental = np.array([4, 7.1, 10.3, 14.9, 40])

    D_NMC622_array = np.linspace(D_NMC622_experimental[0], D_NMC622_experimental[-1], 5000)

    shape_Q_0_NMC622_guess = 10
    loc_Q_0_NMC622_guess = 0
    scale_Q_0_NMC622_guess = 1

    # Q_0_NMC622_guess = f_gamma_cdf(D_NMC622_array, shape_Q_0_NMC622_guess, scale_Q_0_NMC622_guess)
    Q_0_NMC622_guess = stats.gamma.cdf(D_NMC622_array, shape_Q_0_NMC622_guess,loc_Q_0_NMC622_guess, scale_Q_0_NMC622_guess)

    [shape_Q_0_NMC622_opt, loc_Q_0_NMC622_opt,scale_Q_0_NMC622_opt], pcov = optimize.curve_fit(
        f=stats.gamma.cdf,  # model function
        xdata=D_NMC622_experimental,  # x data
        ydata=Q_0_NMC622_experimental,  # y data
        p0=(shape_Q_0_NMC622_guess, loc_Q_0_NMC622_guess, scale_Q_0_NMC622_guess),  # initial value of the parameters
    )
    Q_0_NMC622 = stats.gamma.cdf(D_NMC622_array, shape_Q_0_NMC622_opt, loc_Q_0_NMC622_opt, scale_Q_0_NMC622_opt)
    q_0_NMC622 = stats.gamma.pdf(D_NMC622_array, shape_Q_0_NMC622_opt, loc_Q_0_NMC622_opt, scale_Q_0_NMC622_opt)
    # print(popt)
    print('Gamma parameters for Q_0_NMC622\n' + 'shape_Q_0_NMC622 = ' + str(
        shape_Q_0_NMC622_opt) + '\nloc_Q_0_NMC622 = '
          + str(loc_Q_0_NMC622_opt)+ '\nscale_Q_0_NMC622 = '
          + str(scale_Q_0_NMC622_opt))

    fig_Q_0_NMC622, ax_Q_0_NMC622 = plt.subplots()
    lns_Q_0_NMC622_experimental = ax_Q_0_NMC622.plot(D_NMC622_experimental, Q_0_NMC622_experimental, marker='X',
                                                     markersize=16, linewidth=0, label='Experimental points')
    lns_Q_0_NMC622_guess = ax_Q_0_NMC622.plot(D_NMC622_array, Q_0_NMC622_guess, label=r'$Q_0$ guess')
    lns_Q_0_NMC622_fit = ax_Q_0_NMC622.plot(D_NMC622_array, Q_0_NMC622, label=r'$Q_0$ fit')
    # ax_Q_0_NMC622.set_title('Cumulative distribution Q0 of NMC622')
    ax_Q_0_NMC622.legend(loc='best')
    ax_Q_0_NMC622.set_xlabel('Particle diameter [µm]')
    ax_Q_0_NMC622.set_ylabel('Cumulative density function $Q_0$ $[-]$')

    fig_q_0_NMC622, ax_q_0_NMC622 = plt.subplots()
    # lns_q_3_NMC622_experimental = ax_Q_3_NMC622.plot(D_NMC622_experimental, Q_3_NMC622_experimental, marker='X',
    #                                                  markersize=16, linewidth=0, label='Experimental points')
    # lns_q_3_NMC622_guess = ax_q_3_NMC622.plot(D_NMC622_array, q_3_NMC622_guess, label=r'$Q_3$ guess')
    lns_q_0_NMC622_fit = ax_q_0_NMC622.plot(D_NMC622_array, q_0_NMC622, label=r'$q_0$')
    # ax_q_0_NMC622.set_title('Probability distribution function q0 of NMC622')
    # ax_q_0_NMC622.legend(loc='best')
    ax_q_0_NMC622.set_xlabel('Particle diameter [µm]')
    ax_q_0_NMC622.set_ylabel('Probability density function $q_0$ $[µm^{-1}]$')



    fig_all_distributions, ax_all_distributions = plt.subplots()
    lns_q3_laser_diffraction = ax_all_distributions.plot(10**particle_diameter_exponent_experimental_laser_diffraction,
                                                  q3_lin_experimental_laser_diffraction, label='Laser diffraction - $q_3$')
    lns_q3_tomography = ax_all_distributions.plot(10**particle_diameter_exponent_experimental_tomography,
                                           q_3_lin_experimental_tomography, label=r'Tomography - $q_3$')
    lns_q0_ebner = ax_all_distributions.plot(Ebner_size_array,
                                   Ebner_q0_pdf, label=r'Tomography, raw data - fitted $q_0$')
    lns_q_0_NMC622 = ax_all_distributions.plot(D_NMC622_array, q_0_NMC622, label=r'Purchased NMC-633 $q_0$')
    ax_all_distributions.set_xlabel('particle diameter [µm]')
    ax_all_distributions.set_ylabel('Density distribution $q_i$ $[µm^{-1}]$')
    ax_all_distributions.legend(loc='best')




    fig_data_hist,ax_data_hist = plt.subplots()
    lns_data_hist = ax_data_hist.hist(Ebner_q_0_sample,bins=100,density=True)
    ax_data_hist.set_title('Histogram of created sample')


    D_min = 8
    D_max = 24
    radius_data_for_simulation = 1E-2 * Ebner_q_0_sample[(Ebner_q_0_sample > D_min) & (Ebner_q_0_sample < D_max)] / 2

    # V_part_reduced_radius = np.sum(
    #     4 / 3 * 3.14 * ((1E2 * radius_data_for_simulation[radius_data_for_simulation < (12.073 / 2)]) ** 3))
    #
    # print(V_part_reduced_radius / V_tot)

    print('Max radius in simulation: '+str(max(radius_data_for_simulation)))
    print('Min radius in simulation: ' + str(min(radius_data_for_simulation)))
    np.savetxt("G:/My Drive/Skola/KTH/PhD/Litteratur/X-Ray Tomography of Porous, Transition Metal Oxide Based Lithium Ion Battery Electrodes/Particle_radius_data_from_Ebner_raw_data_gamma_k_284_loc_0457_theta_562_r_min_4_r_max_12.dat",radius_data_for_simulation,fmt='%.8f')

    Previoud_raidus_used_pd = pd.read_csv('G:/My Drive/Skola/KTH/PhD/Litteratur/X-Ray Tomography of Porous, Transition Metal Oxide Based Lithium Ion Battery Electrodes/particle_radius_distr_gamma_k_492_theta_263_r_min_4_r_max_12_cm_size.dat',names=['radius'])
    print(Previoud_raidus_used_pd)
    Previoud_raidus_used = Previoud_raidus_used_pd['radius'].to_numpy()


    fig_input_radius_hist,ax_input_radius_hist = plt.subplots()
    lns_input_radius_hist = ax_input_radius_hist.hist([radius_data_for_simulation,Previoud_raidus_used],label=['new data','old data'], density=True ,bins=100)
    ax_input_radius_hist.set_title('Histogram of particle radius for model input')
    ax_input_radius_hist.legend(loc='best')
    # Eq. 25 in Stoyan and Unland 2022
    # print(D_NMC622_array)
    # print(q_3_NMC622)
    # Q_0_NMC622 = np.array([])
    # for x in D_NMC622_array[1:]:
    #     Q_0_NMC622_integration = integrate.cumtrapz(
    #         D_NMC622_array[D_NMC622_array <= x] ** -3 * D_NMC622_array[D_NMC622_array <= x] ** (shape_Q_3_NMC622_opt - 1) *(np.exp(-D_NMC622_array[D_NMC622_array <= x] / scale_Q_3_NMC622_opt) / (sps.gamma(shape_Q_3_NMC622_opt) * scale_Q_3_NMC622_opt ** shape_Q_3_NMC622_opt)),D_NMC622_array[D_NMC622_array <= x])[-1] / integrate.cumtrapz(
    #                                  D_NMC622_array ** -3 * D_NMC622_array ** (shape_Q_3_NMC622_opt - 1) *
    #                                  (np.exp(-D_NMC622_array / scale_Q_3_NMC622_opt) / (sps.gamma(
    #                                      shape_Q_3_NMC622_opt) * scale_Q_3_NMC622_opt ** shape_Q_3_NMC622_opt)),
    #                                  D_NMC622_array)[-1]
    #     Q_0_NMC622 = np.append(Q_0_NMC622, Q_0_NMC622_integration)
    # Q_0_NMC622 = np.append(0, Q_0_NMC622)
    # q_0_NMC622 = np.gradient(Q_0_NMC622, D_NMC622_array)
    # # q_0_NMC622[0] = 0
    #
    # fig_Q_0_NMC622, ax_Q_0_NMC622 = plt.subplots()
    # # lns_Q_3_NMC622_experimental = ax_Q_3_NMC622.plot(D_NMC622_experimental,Q_3_NMC622_experimental,marker='X',markersize=16,linewidth=0,label='Experimental points')
    # # lns_Q_3_NMC622_guess = ax_Q_3_NMC622.plot(D_NMC622_array,Q_3_NMC622_guess,label=r'$Q_3$ guess')
    # lns_Q_0_NMC622_fit = ax_Q_0_NMC622.plot(D_NMC622_array, Q_0_NMC622, label=r'$Q_0$')
    # ax_Q_0_NMC622.set_title('Probability density function $Q_0$ of NMC622')
    # ax_Q_0_NMC622.legend(loc='best')
    # ax_Q_0_NMC622.set_ylim(ymin=0)
    #
    # D_NMC622_full_array = np.linspace(0, D_NMC622_experimental[-1], 5000)
    # q_0_NMC622_full = f_gamma_pdf(D_NMC622_full_array,shape_Q_3_NMC622_opt,scale_Q_3_NMC622_opt)
    #
    # fig_q_0_NMC622, ax_q_0_NMC622 = plt.subplots()
    # # lns_Q_3_NMC622_experimental = ax_Q_3_NMC622.plot(D_NMC622_experimental,Q_3_NMC622_experimental,marker='X',markersize=16,linewidth=0,label='Experimental points')
    # # lns_Q_3_NMC622_guess = ax_Q_3_NMC622.plot(D_NMC622_array,Q_3_NMC622_guess,label=r'$Q_3$ guess')
    # lns_q_0_NMC622_fit = ax_q_0_NMC622.plot(D_NMC622_full_array, q_0_NMC622_full, label=r'$q_0$')
    # ax_q_0_NMC622.set_title('Probability density function $q_0$ of NMC622')
    # ax_q_0_NMC622.legend(loc='best')
    # ax_q_0_NMC622.set_ylim(ymin=0)






    plt.show()

    # =OLD SCRIPT====================================================================================================
    # shape_q3, scale_q3 = 4.9242, 2.6274  # mean=4, std=2*sqrt(2)
    # max_lim = 24
    # min_lim = 8

    # s_q3 = np.random.gamma(shape_q3, scale_q3, 10000)
#    s = s[ (s>=2*2) & (s<=12*2)]
#     fig_hist,ax_hist = plt.subplots()
#     count, bins, ignored = plt.hist(s_q3, 200, density=True)

#    xarr = np.array(np.linspace(0,50,100))
#    print(type(xarr))
#    y = xarr ** (shape - 1) * (np.exp(-xarr / scale) / (sps.gamma(shape) * scale ** shape))
#     y = bins ** (shape_q3 - 1) * (np.exp(-bins / scale_q3) / (sps.gamma(shape_q3) * scale_q3 ** shape_q3))

# y_int = integrate.cumtrapz(y, bins)
# print(y_int)

#     plt.plot(bins, y, linewidth=2, color='r')
#     plt.xlabel('Particle diameter [µm]', fontsize=16, weight='bold')
# #    plt.ylabel('Number of particles [-]')
#     plt.ylabel('Probability density function [-]', fontsize=16, weight='bold')
#     plt.xlim([4, 24])
#    # plt.legend(loc='best')
#     plt.text(20,.1,s=r'$k=4.92 $'+'\n'+r'$\theta=1.99$' , fontsize=16)
#
#     PDF_fig,ax = plt.subplots()
#     ax.plot(bins, y, linewidth=3)
#     # ax.text(19,.09,s='Gamma\n'+r'$k=4.92 $'+'\n'+r'$\theta=2.63$' , fontsize=16)
#     ax.vlines([min_lim,max_lim],ymin=0,ymax=.1,color='orange')
#     ax.xaxis.set_major_locator(MultipleLocator(5))
#     ax.set_xlim(0, 40)
#     ax.set_ylim(ymin=0,ymax=.1)
#     ax.set_xlabel('Particle diameter [µm]', fontsize=16, weight='bold')
#     #    plt.ylabel('Number of particles [-]')
#     ax.set_ylabel('Probability density function [-]', fontsize=16, weight='bold')
#     ax.tick_params(axis='both', which='major', labelsize=16)
# ==================================================================================================================
