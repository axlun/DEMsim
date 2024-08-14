#
# Created by axlun on 24/08/12
#

import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

plt.style.use('axel_style')

if __name__ == '__main__':
    coefficient_file = "./fitting_coefficient_files/contact_fit_coefficients_E_140_SY_2_9.dat"
    # coefficient_file = "./fitting_coefficient_files/contact_fit_coefficientsE_140_SY_3_3.dat"
    coefficients = pd.read_csv(coefficient_file, header=None).to_numpy(dtype='d')

    print(coefficients[0][0])
    normalised_overlap_range = np.linspace(0, .25, 100)

    test = coefficients[0] + coefficients[1]
    print(test)

    normalised_force = coefficients[0][0] * normalised_overlap_range ** coefficients[1][0] + \
        coefficients[2][0] * normalised_overlap_range ** coefficients[3][0]


    plt.plot(normalised_overlap_range, normalised_force)
    plt.xlabel('h/$R_0$')
    plt.ylabel('F/$R_0^2 \sigma_Y$')

    plt.show()