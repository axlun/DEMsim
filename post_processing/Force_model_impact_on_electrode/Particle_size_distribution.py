import glob
import os.path

from math import pi
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.ticker import (MultipleLocator, ScalarFormatter)

matplotlib.style.use('classic')


def data_grabber(data_directory):
    data_points = np.genfromtxt(data_directory, delimiter=', ')
    #print(data_points[:, 0])
    # data_points = [data_points[:, 1], data_points[:, 2]]
    xdata = data_points[:, 3]
    ydata = data_points[:, 7]
    return xdata, ydata



def main():
    simulation_directory = 'c:/Users/Axel/Documents/DEM/results/electrode_calendaring/bt065N3000'
    file_type = "/particles/*.dou"
    files = glob.glob(simulation_directory + file_type)
    latest_file = max(files, key=os.path.getctime)
    #print(latest_file)
    particle_data = np.genfromtxt(latest_file, delimiter=',')
    #print(particle_data)
    z_position,radius = data_grabber(latest_file)

    z_position = [element * 100 for element in z_position]
    diameter = [element * 2*100 for element in radius]

    print(z_position)
    print(diameter)

    #z_position = [1,85, 90, 2, 4, 40, 46, 67, 70,  99, 3, 5, 6]
    #diameter = [3, 18, 20, 4, 5.5, 6, 10, 12, 16, 22, 21, 20, 21]
    n_x_bins = 50
    n_y_bins =50
    xedges = [np.min(diameter)*.05]
    yedges= [0]
    for i in range(n_x_bins):
        xedges.append(((np.max(diameter))*(i+1)/n_x_bins))
    for i in range(n_y_bins):
        yedges.append(np.max(z_position)*(i+1)/n_y_bins)
    # print(xedges)
    print(yedges)

    xedges = 10 ** np.linspace(np.log10(np.min(diameter)), np.log10(np.max(diameter)),n_x_bins+1)
    #xedges = [4,5,6,7,8,10,12,14,17,20,24]
    # print(xedges)

    H, xedges, yedges = np.histogram2d(diameter,z_position, bins=[xedges, yedges])
    # print(xedges)
    # print(H)
    H=H.T
   # print(H)
    fig = plt.figure(figsize=(7, 7))
    ax = fig.add_subplot(111, title='Particle size distribution of RVE in height direction', xlabel='Particle diameter d (cm)',
                         ylabel='Distancle from current collector z (cm)')
    #plt.imshow(H, interpolation='nearest',aspect="auto", origin='lower', extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])
    im=ax.pcolormesh(xedges, yedges, H)
    ax.set_xscale('log')
    ax.xaxis.set_major_formatter(ScalarFormatter())
    ax.xaxis.set_minor_formatter(ScalarFormatter())
    ax.yaxis.set_major_locator(MultipleLocator(10))
    plt.xlim([3, 30])
    plt.ylim([0, 111])
    # Loop over data dimensions and create text annotations.
    # for i in range((n_x_bins)):
    #     for j in range((n_y_bins)):
    #         text = ax.text((xedges[i+1]-xedges[i])/2+xedges[i], (yedges[j+1]-yedges[j])/2+yedges[j], str(H.T[i, j]),
    #                        ha="center", va="center", color="w")

    # fig2 = plt.figure(2,figsize=(7, 7))
    #
    # plt.hist(diameter, bins=xedges)
    # #H_test, xedges = np.histogram(diameter, bins=xedges)
    # #print(H_test)
    # #plt.imshow(H_test)
    # fig2.gca().set_xscale("log")

    #Colour bar
    cbar = fig.colorbar(im)
    cbar.ax.yaxis.set_major_locator(MultipleLocator(1))
    cbar.ax.set_ylabel('Number of particles', rotation=-90, va="bottom")

    plt.show()


if __name__ == '__main__':
    main()