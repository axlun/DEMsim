import numpy as np

if __name__ == '__main__':

    r_min = 4e-2 # [m]
    mass_scaling = 1e0 # [-]
    rho = 4630 # [kg / m3]
    m_min = rho * 4 * 3.1415 / 3 * r_min**3.0 # [kg]
    R_0 = r_min / 2 # [m]
    E = 140e9 # [Pa]
    nu = .3 # [-]
    E_eff = E / (2 * (1 - nu ** 2)) # [Pa]
    k_p = 4 / 3 * E_eff * R_0**.5 # []
    yield_disp_coeff = 8.59e-3
    h_yield = R_0 * yield_disp_coeff #
    k_hertz = 3 / 2 * k_p * h_yield**.5
    delta_t_min_hertz = 2*(mass_scaling * m_min / k_hertz)**.5

    print("The maximum time step size is " +"{0:.2E}".format(delta_t_min_hertz*1e6)+'µs for the Hertz contact')

    print("Maximum time step with safety factor of 20 is: "+ "{0:.2E}".format(delta_t_min_hertz*1e6/20)+'µs for the Hertz contact' )


    # Rigid plastic contact
    H_bar_max = 2.8
    C_2_max = 1.43
    sigma_Y = 4000e6
    k_rigid_plastic = 3.14*C_2_max*H_bar_max*2*R_0*sigma_Y
    delta_t_min_rigid_plastic = 2 * (mass_scaling * m_min / k_rigid_plastic) ** .5
    print("The maximum time step size is " +"{0:.2E}".format(delta_t_min_rigid_plastic*1e6)+'µs for the rigid plastic '
                                                                                            'contact')

    print("Maximum time step with safety factor of 20 is: "+ "{0:.2E}".format(delta_t_min_rigid_plastic*1e6/20)
          +'µs for the rigid plastic contact' )

