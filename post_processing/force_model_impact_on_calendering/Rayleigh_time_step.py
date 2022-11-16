import numpy as np

if __name__ == '__main__':

    r_min = 4e-2 # [m]
    mass_scaling = 1e2 # [-]
    rho = 4630 # [kg / m3]
    m_min = rho * 4 * 3.1415 / 3 * r_min**3.0 # [kg]
    R_0 = r_min / 2 # [m]
    E = 140e9 # [Pa]
    nu = .3 # [-]
    E_eff = E / (2 * (1 - nu ** 2)) # [Pa]
    k_p = 4 / 3 * E_eff * R_0**.5 # []
    yield_disp_coeff = 8.59e-3
    h_yield = R_0 * yield_disp_coeff #
    k = 3 / 2 * k_p * h_yield**.5
    delta_t_min = 2*(mass_scaling * m_min / k)**.5

#    print("{0:.2E}".format(delta_t_min*1e6))

    print("The maximum time step size is " +"{0:.2E}".format(delta_t_min*1e6)+'µs')

    print("Maximum time step with safety factor of 20 is: "+ "{0:.2E}".format(delta_t_min*1e6/20)+'µs' )