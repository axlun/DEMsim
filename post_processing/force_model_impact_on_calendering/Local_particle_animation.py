from mayavi import mlab

from visualization_functions_3d.animation import Animation, BoundingBox
from visualization_functions_3d.plotting_functions import VelocityArrowPlotter
from tvtk.tools import visual


def main():
    print('start')

#==NATUAL PACKING ======================================================================================================
    # simulation_directory = 'c:/Users/Axel/Documents/DEM/Bertil_results/electrode_natural_packing_hertz/SN_hertz_5000p_btr_8_brr_08_dt_1e0_MS_1e0_el_b_new_tang'
    # simulation_directory = 'c:/Users/Axel/Documents/DEM/Bertil_results/electrode_natural_packing_hertz/SN_hertz_5000p_btr_8_brr_08_dt_1e0_MS_1e0_el_b_new_tang_no_PerBC'

    # simulation_directory = 'c:/Users/Axel/Documents/DEM/Bertil_results/electrode_natural_packing_hertz/SN_hertz_5000p_btr_8_brr_08_dt_1e0_MS_1e0_elast_binder'
    # simulation_directory = 'c:/Users/Axel/Documents/DEM/results/electrode_natural_packing_hertz/SN_hertz_200p_btr_8_brr_08_dt_1e0_MS_1e0_elast_binder_new_tang_W_perBC'
    # simulation_directory = 'c:/Users/Axel/Documents/DEM/Bertil_results//electrode_natural_packing_hertz/SN_hertz_2000p_btr_8_brr_08_dt_1e0_MS_1e0_el_b_new_tang_no_PerBC_fix_gate/'
    # simulation_directory = 'c:/Users/Axel/Documents/DEM/Bertil_results/electrode_natural_packing_hertz/SN_hertz_2000p_btr_8_brr_08_dt_1e0_MS_1e0_el_b_new_tang_fix_gate2/'
    # simulation_directory = 'c:/Users/Axel/Documents/DEM/Bertil_results/electrode_natural_packing_hertz/SN_hertz_2000p_btr_8_brr_08_dt_1e0_MS_1e0_el_b_new_tang_fix_gate3/'
    # simulation_directory = 'c:/Users/Axel/Documents/DEM/Bertil_results/electrode_natural_packing_hertz/SN_hertz_500p_btr_8_brr_08_dt_1e0_MS_1e0_el_b_new_tang_fix_gate4_fullosnign3/'
    # simulation_directory = 'c:/Users/Axel/Documents/DEM/Bertil_results/final_runs/SN_run_1/electrode_natural_packing_hertz/'
    # simulation_directory = 'c:/Users/Axel/Documents/DEM/Bertil_results/final_runs/SN_run_1_q_0_crash/electrode_natural_packing_hertz'
    # simulation_directory = "c:/Users/Axel/Documents/DEM/results/electrode_natural_packing_rigid_plastic_particle/SN_0/"
#==CALENDERING==========================================================================================================
    # simulation_directory = 'c:/Users/Axel/Documents/DEM/Bertil_results/electrode_calendering_hertz/SN_ref_run_1_5000p_btr_5_brr_15_comp_time_20_hal_105_dt_1e2_MS_1e4'
    simulation_directory = 'c:/Users/Axel/Documents/DEM/Bertil_results/final_runs/SN_run_1_rigid_plastic_particle/electrode_natural_packing_rigid_plastic_SY_4GPa_crash/'

#==MECHANICAL LOADING===================================================================================================
    # simulation_directory = 'c:/Users/Axel/Documents/DEM/Bertil_results/electrode_mechanical_loading_hertz/SN_ref_run_1_5000p_btr_5_brr_15_dt_5e1_MS_1e2_SR_2e-3_no_new_binder_run_3_tension'
    # simulation_directory = 'c:/Users/Axel/Documents/DEM/Bertil_results/electrode_mechanical_loading_hertz/SN_ref_run_1_5000p_btr_5_brr_15_dt_5e1_MS_1e2_SR_2e-3_tension'

# ==TANGENTIAL CONTACT TESTING==========================================================================================
#     simulation_directory = "C:/Users/Axel/Documents/DEM/results/contact_testing/elastic_plastic_binder_hertz_particle/tangential_force_wall_contact/New_tangential_force_relation1/"



    fig = mlab.figure(size=(1920, 1200), bgcolor=(1., 1., 1.), fgcolor=(0, 0., 0.))
    visual.set_viewer(fig)
    animation = Animation(simulation_directory)
    animation.save_directory = simulation_directory+'/animation/imgs/'
    animation.save_frames = True


    #Initiate the bounding box
    bbox = BoundingBox()
    #Set limits in the x,y and z direction
    bbox.x_max = lambda t: 2.5
    bbox.x_min = lambda t: -2.5
    bbox.y_max = lambda t: 2.5
    bbox.y_min = lambda t: -2.5
    # bbox.z_max = lambda t: 2.5
    #Apply bounding box to surfaces
    animation.bounding_boxes[0] = bbox
    animation.bounding_boxes[1] = bbox
    animation.bounding_boxes[2] = bbox
    animation.bounding_boxes[3] = bbox
    animation.bounding_boxes[4] = bbox
    animation.bounding_boxes[5] = bbox

    # How long should a frame be shown
    animation.delay = 1e-6

    # animation.surfaces_opacities[0] = 1
    # animation.surfaces_opacities[1] = 1
    #What should be plotted
    animation.plot_force_arrow = False
    # animation.n_force = 50
    # animation.plot_velocity_arrow = True
    # animation.n_vel = 50

    animation.plot_periodic_bc = False
    animation.mirror_particles = False
    animation.view_surfaces = True
    animation.nth_timeframe = 1
    animation.start_time = 0#0.0
    # animation.end_time = 0.01
    # animation.end_time = animation.start_time +0.00011

    animation.run()
    mlab.show()

if __name__ == '__main__':
    main()