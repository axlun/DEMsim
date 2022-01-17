import glob
import os.path

from math import pi
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

matplotlib.style.use('classic')


def data_grabber(data_directory):
    data_points = np.genfromtxt(data_directory, delimiter=',	')
    print(data_points[:, 0])
    # data_points = [data_points[:, 1], data_points[:, 2]]
    ydata = data_points[:, 1]
    xdata = data_points[:, 0]
    return xdata, ydata


def main():
    simulation_directory = 'c:/Users/Axel/Documents/DEM/results/electrode_calendaring/bt065N3000/contacts'
    #simulation_directory = 'c:/Users/Axel/Desktop/contact_test'
    file_type = "/*.dou"
    files = glob.glob(simulation_directory + file_type)
    latest_file = max(files, key=os.path.getctime)
    contact_data = np.genfromtxt(latest_file, delimiter=',')
    #print(contact_data)
    particle_array = [0] * int(contact_data[-1, 0] - contact_data[0, 0] + 1)
    legend = "N = "+str(int(contact_data[-1, 0] - contact_data[0, 0] + 1)) + " particles"
    print(legend)
    for contact in contact_data:
        if contact[5] >= 0:
            particle_array[int(contact[0] - contact_data[0, 0])] += 1
            if contact[1] >= contact_data[0, 0]:
                particle_array[int(contact[1] - contact_data[0, 0])] += 1
    #print(particle_array)
    print(sum(particle_array))
#        if contact[19] != 1 and contact[6] != 0:
#            particle_array[int(contact[0]-contact_data[0, 0])] += 1
#            if contact[1]+1 >= contact_data[0, 0]:
#                particle_array[int(contact[1]-contact_data[0, 0])] += 1
#    print(particle_array)

    # Creating histogram
    fig, ax = plt.subplots(figsize=(10, 7))
    ax.hist(particle_array, bins=range(0,25))
    # Labels
    plt.xlabel("Number of particle-particle contacts")
    plt.ylabel("Number of particles")
    plt.legend([legend])
    plt.title('Histogram of particle-particle contacts')
    # Show plot
    plt.show()
if __name__ == '__main__':
    main()
