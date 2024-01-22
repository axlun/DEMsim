from Bertil_functions.Bertil_functions import *
import os

from mayavi import mlab
from tvtk.tools import visual

from visualization_functions_3d.plotting_functions import BoundingBox
from visualization_functions_3d.battery_contact_plotter import BatteryContactPlotter
from visualization_functions_3d import colors
from force_model_impact_on_calendering.Bertil_snapshot import Snapshot


if __name__ == '__main__':


    plot_contacts = False
    velocity_arrow = False
    force_arrow = False
    n_force = 100
    mirror_particles = True
    periodic_BCs = False


    simulation_directory = "C:/Users/Axel/Documents/DEM/Bertil_results/article_2/viscoelastic_testing/SN_4/electrode_natural_packing_el_pl_binder_el_pl_particle/"
    time_step = 10.65

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
    bbox.x_max = lambda t: 25
    bbox.x_min = lambda t: -25
    bbox.y_max = lambda t: 25
    bbox.y_min = lambda t: -25
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