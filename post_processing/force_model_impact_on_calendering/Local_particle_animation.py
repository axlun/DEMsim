from mayavi import mlab

from visualization_functions_3d.animation import Animation, BoundingBox
from visualization_functions_3d.plotting_functions import VelocityArrowPlotter
from tvtk.tools import visual


def main():
    print('start')

    #==NATUAL PACKING ======================================================================================================
    # simulation_directory = 'c:/Users/Axel/Documents/DEM/Bertil_results/article_2/final_runs/particle_contact_model/SN_2_5/electrode_natural_packing_el_pl_binder_el_pl_particle/'

    #==CALENDERING======================================================================================================
    # simulation_directory = 'c:/Users/Axel/Documents/DEM/Bertil_results/article_2/final_runs/particle_contact_model/SN_2_6/electrode_calendering_el_pl_binder_el_pl_particle/'

    #==MECHANICAL LOADING===============================================================================================
    # simulation_directory = 'c:/Users/Axel/Documents/DEM/Bertil_results/article_2/final_runs/ref_sim/SN_3/electrode_mechanical_loading_hertz_compression/'

    #==RELAXATION=======================================================================================================
    # simulation_directory = 'c:/Users/Axel/Documents/DEM/Bertil_results/article_2/viscoelastic_testing/SN_1/electrode_relaxation_el_pl_binder_el_pl_particle_eps_dot_2e_2_tension/'

    # ==TANGENTIAL CONTACT TESTING==========================================================================================
    # simulation_directory = "C:/Users/Axel/Documents/DEM/results/contact_testing/elastic_plastic_binder_hertz_particle/tangential_force_wall_contact/New_tangential_force_relation1/"

    # ==SWELLING PARTICLE TESTING=======================================================================================
    # simulation_directory = "C:/Users/Axel/Documents/DEM/results/particle_tests/swelling/swelling_cube_die_compaction/no_restart/"

    # ==PERIODIC BC RESTART TEST=======================================================================================
    # simulation_directory = "C:/Users/Axel/Documents/DEM/results/periodic_bc_tests/no_restart/"
    # simulation_directory = "C:/Users/Axel/Documents/DEM/results/periodic_bc_tests/restart/"

    # ==SINTERING=======================================================================================================
    # simulation_directory = "C:/Users/Axel/Documents/DEM/results/sintering/SN_1/"
    # simulation_directory = "C:/Users/Axel/Documents/DEM/results/sintering/rho_50/"
    # simulation_directory = "C:/Users/Axel/Documents/DEM/results/sintering/rho_60/"
    # simulation_directory = "C:/Users/Axel/Documents/DEM/results/sintering/rho_64_5/"

    # ==PERIODIC_BC_TEST================================================================================================
    # simulation_directory = 'c:/Users/Axel/Documents/DEM/results/periodic_bc_tests/no_restart/'

    # ==PERIODIC_PACKING================================================================================================
    simulation_directory = 'c:/Users/Axel/Documents/DEM/results/swelling_electrode/SN_5/electrode_swelling/'

    fig = mlab.figure(size=(1920, 1200), bgcolor=(1., 1., 1.), fgcolor=(0, 0., 0.))
    visual.set_viewer(fig)
    animation = Animation(simulation_directory)
    animation.save_directory = simulation_directory+'/animation/imgs/'
    animation.save_frames = True


    #Initiate the bounding box
    bbox = BoundingBox()
    #Set limits in the x,y and z direction
    # bbox.x_max = lambda t: 250
    # bbox.x_min = lambda t: -250
    # bbox.y_max = lambda t: 250
    # bbox.y_min = lambda t: -250
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

    # animation.surfaces_opacities[0] = .5
    # animation.surfaces_opacities[1] = 1
    # animation.surfaces_opacities[2] = 0
    # animation.surfaces_opacities[3] = 0
    # animation.surfaces_opacities[4] = 0
    # animation.surfaces_opacities[5] = 0

    #What should be plotted
    # animation.plot_force_arrow = False
    # animation.n_force = 50
    # animation.plot_velocity_arrow = True
    # animation.n_vel = 50

    animation.plot_periodic_bc = False
    animation.mirror_particles = True
    animation.view_surfaces = True
    animation.nth_timeframe = 3
    animation.start_time = 0
    # animation.end_time = 0.5

    animation.run()
    mlab.show()

if __name__ == '__main__':
    main()