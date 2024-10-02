from Bertil_functions.Bertil_functions import *
from force_model_impact_on_calendering.Bertil_mechanical_properties import stress_and_linear_strain_finder
import matplotlib.pyplot as plt
import shutil
import os
# from force_model_impact_on_calendering.Bertil_calendering_pressure_multiple_simulations \
#     import calendering_plot_processing
import pandas


class Simulation:
    def __init__(self, sim_dir, plot_contacts = 0):
        self.sim_dir = sim_dir
        self.force_data, self.surface_force_index, self.surface_position_index, \
            self.surface_position_data, self.periodic_BC_data, self.force_fabric_tensor_data, \
            self.ke = bertil_data_gatherer(sim_dir)

        self.time, self.linear_strain, self.sxx, self.syy, self.szz, self.tau_xy, self.tau_xz, self.tau_yx, \
            self.tau_yz, self.tau_zx, self.tau_zy = \
            stress_and_linear_strain_finder(self.periodic_BC_data, self.force_fabric_tensor_data,
                                            self.surface_position_data)
        self.avg_in_plane_compression = - (self.sxx + self.syy) / 2

        self.height_time, self.avg_height = bertil_layer_height(sim_dir)

        if plot_contacts == 1:
            self.contacts_time, self.particle_contact_vec, self.binder_contact_vec, \
                self.binder_particle_contact_vec = contact_counter_bertil(sim_dir)

if __name__ == '__main__':

    contacts = 0
    swelling_material_scaling = Simulation('/scratch/users/axlun/DEMsim/results/article_3/charge_cycling/SN_1/'
                                           'electrode_swelling_material_scaling/', contacts)

    # ==PLOT PARAMETERS=================================================================================================
    fig_dir = 'c:/temp/figures/SOC/'
    try:
        shutil.rmtree(fig_dir)
        os.mkdir(fig_dir)
    except FileNotFoundError:
        print('Directory does not exist')
        os.mkdir(fig_dir)
    except:
        print('Directory could not be removed')
        quit()
    plt.style.use('axel_style')
    # ==================================================================================================================



    # ===FIG 1 CALENDERING SURFACE PRESSURE SURFACE POSITION TIME=======================================================
    fig_SOC, ax_SOC = plt.subplots()
    ax_SOC.set_ylabel("In-plane stress [MPa]")
    ax_SOC.set_xlabel("Time [s]")
    lns_in_plane_stress = ax_SOC.plot(swelling_material_scaling.time,
                                      swelling_material_scaling.avg_in_plane_compression * 1E-6,'C0')
    ax_layer_height = ax_SOC.twinx()
    ax_layer_height.set_ylabel("Active layer thickness [µm]")
    lns_layer_height = ax_layer_height.plot(swelling_material_scaling.time, swelling_material_scaling.avg_height, 'C1')
    # ax_calendering_surface_pressure.set_title('Calendering surface pressure')

    # lns = lns_in_plane_stress + lns_layer_height
    # labs = [l.get_label() for l in lns]
    # ax_SOC.legend(lns, labs, loc='best')
    fig_SOC.tight_layout()

    fname = fig_dir + 'SOC'
    plt.savefig(fname)

    plt.show()
