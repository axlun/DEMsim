import numpy as np

import matplotlib.pyplot as plt
import matplotlib
matplotlib.style.use('classic')


def main():
    #=============================VALUES FOR SMD22============================================
    # strain_in_percent = np.array([-2,-1.25,-1,-0.5,-0.25,0.25,0.5,1,1.25,2])
    # stress_in_mpa = np.array([360,439.5604396,467.6767677,320,277,180,180,200,200,215])
    #
    # strain_in_percent = np.array([-1.25, -1, -0.5, -0.25, 0.25, 0.5, 1, 1.25])
    # stress_in_mpa = np.array([439.5604396, 467.6767677, 320, 277, 180, 180, 200, 200])

    #=============================VALUES FOR Hertz N=3000============================================
    strain_in_percent_Hertz_N3000 = np.array([-2,-1.25,-1,-0.5,-0.25,0.25,0.5,1,1.25,2])
    stress_in_mpa_Hertz_N3000 = np.array([300,280,230,210,170,120,120,130,140,150])

    #=============================VALUES FOR Hertz N=5000============================================
    strain_in_percent_Hertz_N5000 = np.array([-2,-1.25,-1,-0.5,-0.25,0.25,0.5,1,1.25,2])
    stress_in_mpa_Hertz_N5000 = np.array([310,340,390,230,210,130,140,150,160,170])

    #=============================VALUES FOR Hertz N=5000, bsc=2 ============================================
    strain_in_percent_Hertz_N5000_bsc_2 = np.array([-2,-1.25,-1,-0.5,-0.25,0.25,0.5,1,1.25,2])
    stress_in_mpa_Hertz_N5000_bsc_2 = np.array([390,430,260,270,210,140,150,160,160,180])

    #=============================VALUES FOR Hertz N=5000, bt=03 ============================================
    strain_in_percent_Hertz_N5000_bt_03 = np.array([-2,-1.25,-1,-0.5,-0.25,0.25,0.5,1,1.25,2])
    stress_in_mpa_Hertz_N5000_bt_03 = np.array([410,150,180,180,120,50,50,60,60,80])

    #=============================VALUES FOR Hertz N=5000, bt=04 ============================================
    strain_in_percent_Hertz_N5000_bt_04 = np.array([-2,-1.25,-1,-0.5,-0.25,0.25,0.5,1,1.25,2])
    stress_in_mpa_Hertz_N5000_bt_04 = np.array([320,190,190,150,130,50,50,60,70,70])



    fig, ax1 = plt.subplots()
    ax1.set_xlabel("Strain [%]")
    ax1.set_ylabel("E [MPa]")
#    lns1 = ax1.plot(strain_in_percent, stress_in_mpa,'rx-', label='b_t=0.75')
    ax1.set_ylim([0, 500])
    ax1.set_xlim([-2.2, 2.2])

    lns1 = ax1.plot(strain_in_percent_Hertz_N5000_bt_03, stress_in_mpa_Hertz_N5000_bt_03, color = 'red', linestyle = 'dashed', marker = 'x',
    markerfacecolor = 'blue', markersize = 12,markeredgewidth = 3, linewidth=3,label='bt = 0.3')
    lns1 = ax1.plot(strain_in_percent_Hertz_N5000_bt_04, stress_in_mpa_Hertz_N5000_bt_04, color = 'green', linestyle = 'dashed', marker = 'x',
    markerfacecolor = 'blue', markersize = 12,markeredgewidth = 3, linewidth=3,label='bt = 0.4')
    lns2 = ax1.plot(strain_in_percent_Hertz_N5000, stress_in_mpa_Hertz_N5000, color = 'blue', linestyle = 'dashed', marker = 'x',
    markerfacecolor = 'blue', markersize = 12,markeredgewidth = 3, linewidth=3,label='bt = 0.65')
    plt.legend(loc='best')

    plt.show()

if __name__ == '__main__':
    main()