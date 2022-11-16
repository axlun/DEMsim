import os

from mayavi import mlab
from tvtk.tools import visual

from visualization_functions_3d.plotting_functions import BoundingBox
from visualization_functions_3d.battery_contact_plotter import BatteryContactPlotter
from visualization_functions_3d import colors
from force_model_impact_on_calendering.Bertil_snapshot import Snapshot


if __name__ == '__main__':

    time_step =12.57#88.4191#27.0179

    plot_contacts = False
    velocity_arrow = True
    force_arrow = False
    n_force = 100
    mirror_particles = False
    periodic_BCs = False


    # simulation_directory = 'c:/Users/Axel/Documents/DEM/Bertil_results/electrode_mechanical_loading_hertz/SN_hertz_5000p_btr_8_brr_08_large_part_dt_5e1_MS_1e4_SR_1e-3_compression'
    # simulation_directory =  "C:/Users/Axel/Documents/DEM/results/contact_testing/elastic_plastic_binder_hertz_particle/tangential_force/New_tangential_force_relation/"
    # simulation_directory = 'c:/Users/Axel/Documents/DEM/Bertil_results/electrode_calendering_hertz/SN_hertz_5000p_btr_8_brr_08_comp_time_40_hal_110_dt_1e2_MS_1e4_elast_binder/'
    simulation_directory = "c:/Users/Axel/Documents/DEM/Bertil_results/electrode_natural_packing_hertz/SN_hertz_5000p_btr_8_brr_08_dt_1e0_MS_1e0/"


    fig = mlab.figure(size=(1920, 1200), bgcolor=(1., 1., 1.), fgcolor=(0, 0., 0.))
    visual.set_viewer(fig)

    if plot_contacts:
        snapshot = Snapshot(os.path.expanduser(simulation_directory), BatteryContactPlotter)
        snapshot.contact_plotter.color = colors.white
        snapshot.contact_plotter.binder_radius = 0.02
    else:
        snapshot = Snapshot(os.path.expanduser(simulation_directory))

    snapshot.particle_opacity = 1
    snapshot.mirror_particles_plotter.opacity = 1
    snapshot.mirror_particles = mirror_particles
    snapshot.plot_periodic_bc=periodic_BCs
    snapshot.plot_velocity_arrow = velocity_arrow
    snapshot.plot_force_arrow=force_arrow
    snapshot.n_force = n_force

    # snapshot.surfaces_colors[1] = colors.red
    # snapshot.surfaces_opacities[1] = 0
    # snapshot.surfaces_opacities[0] = 0
    # Initiate the bounding box
    bbox = BoundingBox()
    # Set limits in the x,y and z direction
    bbox.x_max = lambda t: 2.5
    bbox.x_min = lambda t: -2.5
    bbox.y_max = lambda t: 2.5
    bbox.y_min = lambda t: -2.5
    # bbox.z_max = lambda t: 2
    # Apply bounding box to surfaces
    snapshot.surface_bounding_boxes[0] = bbox
    snapshot.surface_bounding_boxes[1] = bbox
    snapshot.surface_bounding_boxes[2] = bbox
    snapshot.surface_bounding_boxes[3] = bbox
    snapshot.surface_bounding_boxes[4] = bbox
    snapshot.surface_bounding_boxes[5] = bbox
    snapshot.plot(time_step)

    print('showing plot')
    mlab.show()