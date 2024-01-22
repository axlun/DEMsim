
import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
import matplotlib
np.set_printoptions(threshold=np.inf)

def force_overlap_fun(R_1,R_2):
    C_1 = 0.7690286505871559
    C_2 = 0.5567414886513075
    C_3 = 0.33897474117487203
    C_4 = 0.15183763225672245

    c_2_max = .84
    H_bar_max = 2.44

    # C_1, C_2, C_3, C_4 = [0.6487628707218379, 0.6004860105394244, 0.36316268485256603, 0.14363844104018408]
    E = 140e9
    nu = .3
    sigma_y = 4e9

    # C_1 = 1.07
    # C_2 = 0.51
    # C_3 = 0.16
    # C_4 = 0.20
    # # C_1, C_2, C_3, C_4 = [0.6487628707218379, 0.6004860105394244, 0.36316268485256603, 0.14363844104018408]
    # E = 200e9
    # nu = .3
    # sigma_y = .2e9


    E_eff = E/(1-nu**2)
    R_eff = (R_1*R_2)/(R_1+R_2)
    R_01 = 2*R_1*R_2/(R_1+R_2)
    # a = np.linspace(sigma_y*R_01/E_eff*np.exp(-C_3/C_4), .5*R_01, 10000)
    print(3*R_01*sigma_y/(E_eff))
    a = np.linspace(3*R_01*sigma_y/E_eff, .8 * R_01, 1000)
    # for x in a:
    #     if (C_3+C_4*np.log(E_eff/sigma_y*x/R_01)) > c_2_max:
    #         print('c^2_max')

    overlap_vec = np.where((C_3+C_4*np.log(E_eff/sigma_y*a/R_01)) < c_2_max, a**2 / (2*R_01) * 1/(C_3+C_4*np.log(E_eff/sigma_y*a/R_01)), a**2 / (2*R_01) * 1/(c_2_max))
    # overlap_vec = a**2 / (2*R_01) * 1/(C_3+C_4*np.log(E_eff/sigma_y*a/R_01))
    overlap_vec_norm = overlap_vec / R_eff

    # for x in a:
    #     if (C_1 + C_2 * np.log(E_eff / sigma_y * x / R_01))>H_bar_max:
    #         print('H_bar_max')

    force_vec = np.where((C_1 + C_2 * np.log(E_eff / sigma_y * a / R_01))<H_bar_max,(C_1+C_2*np.log(E_eff/sigma_y * a/R_01))*sigma_y*a**2,(H_bar_max)*sigma_y*a**2)
    # force_vec = (C_1+C_2*np.log(E_eff/sigma_y * a/R_01))*sigma_y*a**2
    force_vec_norm = force_vec / (R_eff**2 * sigma_y)

    # fig_h_a, ax_h_a = plt.subplots()
    # lns_h_a = ax_h_a.plot(a, overlap_vec_norm,label='overlap norm')
    # lns_F_a = ax_h_a.plot(a,force_vec_norm,label='force norm')
    # ax_h_a.legend(loc='best')
    #
    # fig_F_h, ax_F_h = plt.subplots()
    # lns_F_h = ax_F_h.plot(overlap_vec_norm,force_vec_norm,label='Force-overlap')
    # # lns_F_a = ax_h_a.plot(a,force_vec_norm,label='force norm')
    # ax_F_h.legend(loc='best')
    return overlap_vec, force_vec, overlap_vec_norm, force_vec_norm


def rigid_plast_force_overlap_fun(R_1,R_2):
    E = 140e9
    nu = .3
    sigma_y = 4e9
    sigma_y_eff = 1/(1/sigma_y + 1/sigma_y)

    c_2_max = .84
    H_bar_max = 2.45
    # E = 200e9
    # nu = .3
    # sigma_y = .2e9

    E_eff = E /2/ (1 - nu ** 2)
    R_eff = (R_1 * R_2) / (R_1 + R_2)
    overlap_vec = np.linspace(0,0.5*R_eff,1000)
    overlap_vec_norm = overlap_vec /R_eff
    force_vec = 2*3*3.14*c_2_max * sigma_y_eff * R_eff * overlap_vec
    force_vec_norm = force_vec/(R_eff**2 * sigma_y)
    return overlap_vec, force_vec, overlap_vec_norm, force_vec_norm


def hertz_force_overlap_fun(R_1,R_2):
    E = 140e9
    nu = .3
    sigma_y = 4e9

    # E = 200e9
    # nu = .3
    # sigma_y = .2e9

    E_eff = E /2/ (1 - nu ** 2)
    R_eff = (R_1 * R_2) / (R_1 + R_2)
    overlap_vec = np.linspace(0,0.5*R_eff,1000)
    overlap_vec_norm = overlap_vec /R_eff
    force_vec = 4/3 * E_eff * R_eff**(1/2) * overlap_vec**(3/2)
    force_vec_norm = force_vec/(R_eff**2 * sigma_y)
    return overlap_vec, force_vec, overlap_vec_norm, force_vec_norm

def hertz_plast_force_overlap_fun(R_1,R_2):
    E = 140e9
    nu = .3
    sigma_y = 4e9
    delta_y = 8.59e-3

    # E = 200e9
    # nu = .3
    # sigma_y = .2e9

    E_eff = E /2/ (1 - nu ** 2)
    R_eff = (R_1 * R_2) / (R_1 + R_2)
    overlap_vec = np.linspace(0,0.5*R_eff,1000)
    overlap_vec_norm = overlap_vec /R_eff
    force_vec = np.array([])
    for x in overlap_vec:
        if x > delta_y*R_eff:
            force_vec = np.append(force_vec, 4/3 * E_eff * R_eff**(1/2) * (delta_y*R_eff)**(3/2) + 2*E_eff*(R_eff*delta_y*R_eff)**.5*(x-delta_y*R_eff))
        else:
            force_vec = np.append(force_vec,4/3 * E_eff * R_eff**(1/2) * x**(3/2))
    # force_vec = 4/3 * E_eff * R_eff**(1/2) * overlap_vec**(3/2)
    force_vec_norm = force_vec/(R_eff**2 * sigma_y)
    return overlap_vec, force_vec, overlap_vec_norm, force_vec_norm


if __name__ == '__main__':
    plt.style.use('axel_style')

    R_mean = 5e-6
    R_min = 4e-6
    R_max = 12e-6

    R_unit = 1

    overlap_vec_mean_mean, force_vec_mean_mean, overlap_vec_norm_mean_mean, force_vec_norm_mean_mean = \
        force_overlap_fun(R_mean, R_mean)
    hertz_overlap_vec_mean_mean, hertz_force_vec_mean_mean, hertz_overlap_vec_norm_mean_mean, hertz_force_vec_norm_mean_mean = \
        hertz_force_overlap_fun(R_mean, R_mean)
    hertz_plast_overlap_vec_mean_mean, hertz_plast_force_vec_mean_mean, hertz_plast_overlap_vec_norm_mean_mean, hertz_plast_force_vec_norm_mean_mean = \
        hertz_plast_force_overlap_fun(R_mean, R_mean)
    rigid_plast_overlap_vec_mean_mean, rigid_plast_force_vec_mean_mean, rigid_plast_overlap_vec_norm_mean_mean, rigid_plast_force_vec_norm_mean_mean = \
        rigid_plast_force_overlap_fun(R_mean, R_mean)

    overlap_vec_max_max, force_vec_max_max, overlap_vec_norm_max_max, force_vec_norm_max_max = \
        force_overlap_fun(R_max, R_max)
    hertz_overlap_vec_max_max, hertz_force_vec_max_max, hertz_overlap_vec_norm_max_max, hertz_force_vec_norm_max_max = \
        hertz_force_overlap_fun(R_max, R_max)

    overlap_vec_min_min, force_vec_min_min, overlap_vec_norm_min_min, force_vec_norm_min_min = \
        force_overlap_fun(R_min, R_min)
    hertz_overlap_vec_min_min, hertz_force_vec_min_min, hertz_overlap_vec_norm_min_min, hertz_force_vec_norm_min_min = \
        hertz_force_overlap_fun(R_min, R_min)

    overlap_vec_min_max, force_vec_min_max, overlap_vec_norm_min_max, force_vec_norm_min_max = \
        force_overlap_fun(R_max, R_min)
    hertz_overlap_vec_min_max, hertz_force_vec_min_max, hertz_overlap_vec_norm_min_max, hertz_force_vec_norm_min_max = \
        hertz_force_overlap_fun(R_max, R_min)

    overlap_vec_unit, force_vec_unit, overlap_vec_norm_unit, force_vec_norm_unit = \
        force_overlap_fun(R_unit, R_unit)
    hertz_overlap_vec_unit, hertz_force_vec_unit, hertz_overlap_vec_norm_unit, hertz_force_vec_norm_unit = \
        hertz_force_overlap_fun(R_unit, R_unit)


    hertz_end_index = np.where(hertz_overlap_vec_norm_unit >= overlap_vec_norm_unit[0])[0][0]-1
    print(hertz_end_index)
    master_overlap_vec_unit = np.append(hertz_overlap_vec_unit[:hertz_end_index],overlap_vec_unit)
    master_force_vec_unit = np.append(hertz_force_vec_unit[:hertz_end_index],force_vec_unit)
    master_overlap_vec_norm_unit = np.append(hertz_overlap_vec_norm_unit[:hertz_end_index],overlap_vec_norm_unit)
    master_force_vec_norm_unit = np.append(hertz_force_vec_norm_unit[:hertz_end_index],force_vec_norm_unit)

    master_force_grad_vec_norm_unit = np.gradient(master_force_vec_norm_unit,master_overlap_vec_norm_unit).T



    # =INTERPOLATION OF FORCE OVERLAP CURVE=============================================================================
    master_norm_force_overlap_fun = interpolate.interp1d(master_overlap_vec_norm_unit,master_force_vec_norm_unit, kind='quadratic', fill_value='extrapolate')
    master_norm_overlap = np.linspace(0,2,10000)
    master_norm_force = master_norm_force_overlap_fun(master_norm_overlap)
    # ==================================================================================================================
    # =INTERPOLATION OF FORCE GRADIENT OVERLAP CURVE====================================================================
    master_norm_force_grad_overlap_fun = interpolate.interp1d(master_overlap_vec_norm_unit, master_force_grad_vec_norm_unit,
                                                         kind='quadratic', fill_value='extrapolate')
    master_norm_force_grad = master_norm_force_grad_overlap_fun(master_norm_overlap)
    # ==================================================================================================================

    fig,ax = plt.subplots()
    lns_hertz_mean_mean = ax.plot(hertz_overlap_vec_norm_mean_mean, hertz_force_vec_norm_mean_mean, label='hertz_mean_mean')
    lns_mean_mean = ax.plot(overlap_vec_norm_mean_mean, force_vec_norm_mean_mean,label ='mean_mean')

    # lns_hertz_plast_mean_mean = ax.plot(hertz_plast_overlap_vec_norm_mean_mean, hertz_plast_force_vec_norm_mean_mean, label='hertz_plast_mean_mean')
    # lns_rigid_plast_mean_mean = ax.plot(rigid_plast_overlap_vec_norm_mean_mean, rigid_plast_force_vec_norm_mean_mean, label='rigid_plast_mean_mean')
    lns_hertz_max_max = ax.plot(hertz_overlap_vec_norm_max_max, hertz_force_vec_norm_max_max,
                                  label='hertz_max_max')
    lns_max_max = ax.plot(overlap_vec_norm_max_max, force_vec_norm_max_max,label ='max_max')

    lns_hertz_min_min = ax.plot(hertz_overlap_vec_norm_min_min, hertz_force_vec_norm_min_min,
                                label='hertz_min_min')
    lns_min_min = ax.plot(overlap_vec_norm_min_min, force_vec_norm_min_min,label ='min_min')

    lns_hertz_min_max = ax.plot(hertz_overlap_vec_norm_min_max, hertz_force_vec_norm_min_max,
                                label='hertz_min_max')
    lns_min_max = ax.plot(overlap_vec_norm_min_max, force_vec_norm_min_max,label ='min_max')

    lns_hertz_unit = ax.plot(hertz_overlap_vec_norm_unit, hertz_force_vec_norm_unit,
                                label='hertz_unit')
    lns_unit = ax.plot(overlap_vec_norm_unit, force_vec_norm_unit,label='unit')

    lns_master = ax.plot(master_overlap_vec_norm_unit,master_force_vec_norm_unit,label = 'master')

    lns_master_interpol = ax.plot(master_norm_overlap,master_norm_force, label='master_interpol')

    ax.set_xlabel(r'$h/R^*$')
    ax.set_ylabel(r'$F/{R^*} ^2 \sigma_Y$')
    ax.legend(loc='best')

    fig_master_curve, ax_master_curve = plt.subplots()
    lns_master = ax_master_curve.plot(master_overlap_vec_norm_unit,master_force_vec_norm_unit,'-*',label = 'master')
    lns_master_interpol = ax_master_curve.plot(master_norm_overlap,master_norm_force, label='master extrapolated')
    lns_hertz_lim = ax_master_curve.plot([hertz_overlap_vec_norm_unit[hertz_end_index-1],hertz_overlap_vec_norm_unit[hertz_end_index-1]],[0,1.2*max(master_norm_force)],label='Hertz limit')
    ax_master_curve.set_xlim(xmin=0)
    ax_master_curve.set_ylim(ymin=0,ymax=1.1*max(master_norm_force))
    ax_master_curve.set_xlabel(r'$h/R^*$')
    ax_master_curve.set_ylabel(r'$F/{R^*} ^2 \sigma_Y$')
    ax_master_curve.legend(loc='best')

    fig_master_curve_slope, ax_master_curve_slope = plt.subplots()
    lns_master_slope = ax_master_curve_slope.plot(master_overlap_vec_norm_unit,master_force_grad_vec_norm_unit,'-*',label = 'master')
    lns_master_interpol_slope = ax_master_curve_slope.plot(master_norm_overlap,master_norm_force_grad, label='master extrapolated')
    lns_hertz_lim_slope = ax_master_curve_slope.plot([hertz_overlap_vec_norm_unit[hertz_end_index-1],hertz_overlap_vec_norm_unit[hertz_end_index-1]],[0,1.2*max(master_norm_force)],label='Hertz limit')
    ax_master_curve_slope.set_xlim(xmin=0)
    ax_master_curve_slope.set_ylim(ymin=0,ymax=1.1*max(master_norm_force_grad))
    ax_master_curve_slope.set_xlabel(r'$h/R^*$')
    ax_master_curve_slope.set_ylabel(r'$F\'/{R^*} ^2 \sigma_Y$')
    ax_master_curve_slope.legend(loc='best')


    plt.show()