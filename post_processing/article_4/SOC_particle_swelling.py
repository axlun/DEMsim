import matplotlib
import matplotlib.pyplot as plt
import numpy as np

plt.style.use('axel_style')

if __name__ == '__main__':
    print('start script')

    x_array = np.arange(0, 7, 1)
    print(x_array)
    soc = np.append(np.tile([0, 100], 3), 0)
    print(soc)
    volume_ratio = np.array(
        [1, 1.2, 1.05, 1.2, 1.066666667, 1.2, 1.083333333])
    print(volume_ratio)
    # material_scaling = np.array(
    #     [1, 0.8, 0.9, 0.8, 0.896551724, 0.8, 0.893103448, 0.8, 0.889655172, 0.8, 0.886206897, 0.8, 0.882758621, 0.8,
    #      0.879310345, 0.8, 0.875862069, 0.8, 0.872413793, 0.8, 0.868965517])
    fig_table, ax_table = plt.subplots()
    lns_vol_ratio = ax_table.plot(x_array, volume_ratio, color='C0', label='Volume ratio')
    ax_table.set_xlim(xmin=0, xmax=3)
    ax_table.set_xlabel('SOC [%]')
    ax_table.set_ylabel('Volume ratio $V/V_0$ [-]')
    # ax_table.set_xlim(xmin=0, xmax=100)
    ymin = 0.95
    ymax = 1.25
    ax_table.set_ylim(ymin=ymin, ymax=ymax)
    tick_array = np.linspace(0, 6, 7)
    tick_label_array = [0, 100, 0, 100, 0, 100, 0]
    ax_table.set_xticks(tick_array)
    ax_table.set_xticklabels(tick_label_array)

    # ax_material_scaling = ax_table.twinx()
    # ax_material_scaling.set_ylabel('Material scaling $k$ [-]')
    # ax_material_scaling.set_ylim(ymin=ymin, ymax=ymax)
    # ax_material_scaling.plot(x_array, material_scaling)
    # lns_mat_scale = ax_material_scaling.plot(x_array, volume_ratio, color='C1', label='Material scaling')
    # lns_mat_scale = ax_material_scaling.plot(x_array, material_scaling, color='C1', label='Material scaling')

    ax_2nd_x = ax_table.twiny()
    ax_2nd_x.plot([0, 6], [0, 0])
    ax_2nd_x.set_xticks(np.linspace(0, 3, 4))
    ax_2nd_x.set_xticklabels(np.arange(0, 4, 1))
    ax_2nd_x.set_xlim(xmin=0, xmax=3)
    ax_2nd_x.set_xlabel('Charge cycle [-]')
    # lns =  lns_vol_ratio + lns_mat_scale
    # labs = [l.get_label() for l in lns]
    # ax_table.legend(lns, labs, loc='upper right')

    fig_table.tight_layout()

    """
    x_array = np.arange(0, 21, 1)
    soc = np.append(np.tile([0, 100], 10), 0)
    print(soc)
    volume_ratio = np.array(
        [1, 1.2, 1.05, 1.2, 1.066666667, 1.2, 1.083333333, 1.2, 1.1, 1.2, 1.116666667, 1.2, 1.133333333, 1.2, 1.15, 1.2,
         1.166666667, 1.2, 1.183333333, 1.2, 1.2])
    print(volume_ratio)
    material_scaling = np.array(
        [1, 0.8, 0.9, 0.8, 0.896551724, 0.8, 0.893103448, 0.8, 0.889655172, 0.8, 0.886206897, 0.8, 0.882758621, 0.8,
         0.879310345, 0.8, 0.875862069, 0.8, 0.872413793, 0.8, 0.868965517])
         
    fig_table, ax_table = plt.subplots()
    lns_vol_ratio = ax_table.plot(x_array, volume_ratio, color='C0', label='Volume ratio')
    ax_table.set_xlim(xmin=0, xmax=20)
    ax_table.set_xlabel('SOC [%]')
    ax_table.set_ylabel('Volume ratio $V/V_0$ [-]')
    # ax_table.set_xlim(xmin=0, xmax=100)
    ymin = 0.75
    ymax = 1.28
    ax_table.set_ylim(ymin=ymin, ymax=ymax)
    tick_array = np.linspace(0, 20, 21)
    tick_label_array = [0, 100, 0, 100, 0, 100, 0, 100, 0, 100, 0, 100, 0, 100, 0, 100, 0, 100, 0, 100, 0]
    ax_table.set_xticks(tick_array)
    ax_table.set_xticklabels(tick_label_array)

    ax_material_scaling = ax_table.twinx()
    ax_material_scaling.set_ylabel('Material scaling $k$ [-]')
    ax_material_scaling.set_ylim(ymin=ymin, ymax=ymax)
    # ax_material_scaling.plot(x_array, material_scaling)
    # lns_mat_scale = ax_material_scaling.plot(x_array, volume_ratio, color='C1', label='Material scaling')
    lns_mat_scale = ax_material_scaling.plot(x_array, material_scaling, color='C1', label='Material scaling')
    ax_2nd_x = ax_table.twiny()
    ax_2nd_x.plot([0, 20], [0, 0])
    ax_2nd_x.set_xticks(np.linspace(0, 20, 11))
    ax_2nd_x.set_xticklabels(np.arange(0, 11, 1))
    ax_2nd_x.set_xlim(xmin=0, xmax=20)
    ax_2nd_x.set_xlabel('Charge cycle [-]')
    lns =  lns_vol_ratio + lns_mat_scale
    labs = [l.get_label() for l in lns]
    ax_table.legend(lns, labs, loc='upper right')

    fig_table.tight_layout()

    # ===================================================================================================================
    fig_table_2, ax_table_2 = plt.subplots(2, 1)
    ax_table_2[1].set_xlabel('SOC [%]')
    ax_table_2[1].set_ylabel('Material degradation factor [-]')

    ax_table_2[0].set_xlabel('Charge cycle [-]')
    ax_table_2[0].set_ylabel('Volume ratio $V/V_0$ [-]')
    ax_table_2[0].xaxis.set_label_position('top')
    ax_table_2[0].plot(x_array, volume_ratio)
    ax_table_2[0].tick_params(top=True, labeltop=True, bottom=False, labelbottom=False)
    tick_array_cycle = np.linspace(0, 20, 11)
    ax_table_2[0].set_xticks(tick_array_cycle)
    tick_label_array_cycle = np.arange(0, 11, 1)
    ax_table_2[0].set_xticklabels(tick_label_array_cycle)
    ax_table_2[1].plot(x_array, material_scaling)
    ax_table_2[1].set_xticks(tick_array)
    ax_table_2[1].set_xticklabels(tick_label_array)

    fig_table_2.tight_layout()
    """
    plt.show()
