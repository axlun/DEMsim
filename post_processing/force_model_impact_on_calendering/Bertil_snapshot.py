from Bertil_functions.Bertil_functions import *
import pathlib
from collections import defaultdict
import os
import shutil
import re

from mayavi import mlab
from tvtk.tools import visual
import numpy as np

from visualization_functions_3d.plotting_functions import SpheresPlotter, SurfacesPlotter, BoundingBox, VelocityArrowPlotter, ForceArrowPlotter
from visualization_functions_3d.periodic_bc import PeriodicBC
from visualization_functions_3d.battery_contact_plotter import BatteryContactPlotter
from visualization_functions_3d import colors


class Snapshot:
    def __init__(self, directory, contact_plotter_class=None):
        directory = pathlib.Path(directory)
        self.directory = directory
        self.plot_periodic_bc = False
        self.periodic_bc_plotter = None
        if os.path.isfile(self.directory / 'periodic_bc.dou'):
            self.periodic_bc_plotter = PeriodicBC(self.directory / 'periodic_bc.dou')
        self.mirror_particles = False
        self.particle_bounding_box = BoundingBox()
        self.surface_bounding_boxes = defaultdict(BoundingBox)
        self.particle_opacity = 1.

        self.surfaces_colors = defaultdict(lambda: colors.blue)
        self.surfaces_opacities = defaultdict(lambda: 0.5)
        self.plot_order = None
        self.visible_functions = defaultdict(lambda: lambda t: True)

        self.spheres_plotter = SpheresPlotter()
        self.plot_velocity_arrow = False
        self.velocity_arrow_plotter = VelocityArrowPlotter()

        self.plot_force_arrow = False
        self.n_force = -1
        # self.force_arrow_plotter = ForceArrowPlotter()


        self.mirror_particles_plotter = SpheresPlotter(color=colors.silver)
        self.surfaces_plotter = None
        if os.path.isfile(self.directory / 'surface_positions.dou'):
            self.surfaces_plotter = SurfacesPlotter(self.directory / 'surface_positions.dou', self.surfaces_colors,
                                                    self.surfaces_opacities, self.plot_order,
                                                    self.surface_bounding_boxes, self.visible_functions)
        if contact_plotter_class is not None:
            self.contact_plotter = contact_plotter_class(directory=str(directory))
        else:
            self.contact_plotter = None

    def create_periodic_bc_plotter(self):
        self.periodic_bc_plotter = PeriodicBC(self.directory / 'periodic_bc.dou')

    def plot(self, time):
        if time.is_integer():
            time = int(time)
        self.spheres_plotter.opacity = self.particle_opacity
        particle_data = np.genfromtxt(self.directory / ('particles/particles_' + str(time) + '.dou'), delimiter=',')
        surface_positions_data = np.genfromtxt(self.directory / ('surface_positions.dou'), delimiter=',')
        self.spheres_plotter.bounding_box = self.particle_bounding_box
        self.spheres_plotter.plot(particle_data)

        if self.plot_velocity_arrow:
            self.velocity_arrow_plotter.plot(self.directory,time,particle_data,surface_positions_data)

        if self.plot_force_arrow:
            self.force_arrow_plotter = ForceArrowPlotter(colors.green, 1., self.n_force)
            self.force_arrow_plotter.plot(particle_data)

        if self.mirror_particles:
            mirror_particle_data = np.genfromtxt(self.directory / ('mirror_particles/mirror_particles_'
                                                                   + str(time) + '.dou'), delimiter=',')
            self.mirror_particles_plotter.plot(mirror_particle_data)
        if self.surfaces_plotter:
            self.surfaces_plotter.bounding_boxes = self.surface_bounding_boxes
            self.surfaces_plotter.surfaces_colors = self.surfaces_colors
            self.surfaces_plotter.surfaces_opacities = self.surfaces_opacities
            self.surfaces_plotter.plot(time)
        if self.plot_periodic_bc:
            self.periodic_bc_plotter.plot(time)
        if self.contact_plotter:
            self.contact_plotter.plot(time)



        f = mlab.gcf()
        f.scene.render()

# NEEDS TO BE FIXED
# GET NEW FILES WHEN SIM IS STILL RUNNING================================================================================
# When getting neighbouring particle files, also get neighbouring mirror_particles etc.
if __name__ == '__main__':

    time_step = 10.9
    velocity_arrow = True
    force_arrow = False
    plot_contacts = False
    mirror_particles = True
    periodic_BCs = False

    # ==NATUAL PACKING =================================================================================================
    simulation_directory = '/scratch/users/axlun/DEMsim/results/electrode_natural_packing_hertz/SN_hertz_5000p_btr_8_brr_05_dt_1e0_MS_1e0/'

    # ==CALENDERING=====================================================================================================
#     simulation_directory = '/scratch/users/axlun/DEMsim/results/electrode_calendering_hertz/SN_hertz_5000p_btr_8_brr_08_comp_time_20_hal_105_dt_1e2_MS_1e4/'


    # ==MECHANICAL LOADING==============================================================================================
#    bertil_simulation_directory = '/scratch/users/axlun/DEMsim/results/electrode_mechanical_loading_hertz/SN_hertz_5000p_btr_8_brr_08_large_part_dt_5e1_MS_1e4_SR_1e-3_compression'

    # ==RESTING=========================================================================================================
#    bertil_simulation_directory = '/scratch/users/axlun/DEMsim/results/electrode_resting_hertz/SN_hertz_5000p_btr_8_brr_08_dt_5e1_MS_1e4_RT_10'

    if (one_file_reader(simulation_directory+'/particles/particles_'+str(time_step)+'.dou') == OSError):
        print("Particle file not existing on Bertil")
        exit()



    # =========================BUILDING TEMPORARY CATALOGE==============================================================

    local_directory = 'c:/Users/Axel/Documents/DEM/Bertil_results/' + re.split('/scratch/users/axlun/DEMsim/results/',
                                                                               simulation_directory)[-1]+'/'
    particles_dir = local_directory + 'particles/'
    mirror_particles_dir = local_directory + 'mirror_particles/'
    contacts_dir = local_directory + 'contacts/'

    if (1==1):# Always update results files #not os.path.exists(particles_dir+ "particles_"+str(time_step)+".dou"):
        print("Result files not on local machine")
        if not os.path.isdir(local_directory):
            print('Result directory not on local machine')
            try:
                os.mkdir(local_directory)
            except OSError:
                print("Creation of the directory %s failed" % local_directory)
            else:
                print("Successfully created the directory %s " % local_directory)
            os.mkdir(particles_dir)
            os.mkdir(contacts_dir)
            os.mkdir(mirror_particles_dir)
        else: print('Result directory is on local machine')

        # check if mirror particles exist

        # =========WRITING SURFACE POSITIONS============================================================================
        if (1==1):#not os.path.exists(local_directory + "surface_positions.dou"):
            print("Writing surface positions file")
            with open(local_directory + "surface_positions.dou", 'w') as fp:
                fp.write(''.join((one_file_reader(simulation_directory + '/surface_positions.dou'))))
        else:
            print('surface_positions  already local')
        # =========WRITING PARTICLES====================================================================================
        print("Writing particle file")

        #also get neighbouring particle time files
        surface_position_data = np.genfromtxt(local_directory + "surface_positions.dou", delimiter=',')
        time_step_index = np.where(surface_position_data[:,-1] == time_step)[0][0]
        prev_time_step_index = time_step_index - 1
        print(prev_time_step_index)
        next_time_step_index = time_step_index + 1
        if prev_time_step_index < 0:
            prev_time_step_index = time_step_index

        if next_time_step_index > np.size(surface_position_data[:,-1])-1:
            next_time_step_index = [np.size(surface_position_data[:,-1])-1]
        #calculateing time increment
        next_time_step = surface_position_data[next_time_step_index,-1]

        if next_time_step.is_integer():
            int(next_time_step)

        previous_time_step = surface_position_data[prev_time_step_index,-1]

        if previous_time_step.is_integer():
            int(previous_time_step)

        print(time_step)
        print(previous_time_step)
        print(next_time_step)
        #Write current particle time
        with open(particles_dir + "particles_"+str(time_step)+".dou", 'w') as fp:
            fp.write(''.join((one_file_reader(simulation_directory+'/particles/particles_'+str(time_step)+'.dou'))))

        #Write previous particle_time
        if not os.path.exists(particles_dir + "particles_" + str(previous_time_step) + ".dou"):
            with open(particles_dir + "particles_"+str(previous_time_step)+".dou", 'w') as fp:
                fp.write(''.join((one_file_reader(simulation_directory+'/particles/particles_'+str(previous_time_step)+'.dou'))))

        # Write next particle_time
        if not os.path.exists(particles_dir + "particles_" + str(next_time_step) + ".dou"):
            with open(particles_dir + "particles_"+str(next_time_step)+".dou", 'w') as fp:
                fp.write(''.join((one_file_reader(simulation_directory+'/particles/particles_'+str(next_time_step)+'.dou'))))

        # =========WRITING CONTACTS=====================================================================================
        if plot_contacts:
            print("Writing contact file")
            with open(contacts_dir + "contacts_"+str(time_step)+".dou", 'w') as fp:
                fp.write(''.join((one_file_reader(simulation_directory+'/contacts/contacts_'+str(time_step)+'.dou'))))

            # Write previous particle_time
            if not os.path.exists(contacts_dir + "contacts_" + str(previous_time_step) + ".dou"):
                with open(contacts_dir + "contacts_" + str(previous_time_step) + ".dou", 'w') as fp:
                    fp.write(''.join((one_file_reader(
                        simulation_directory + '/contacts/contacts_' + str(previous_time_step) + '.dou'))))

            # Write next particle_time
            if not os.path.exists(contacts_dir + "contacts_" + str(next_time_step) + ".dou"):
                with open(contacts_dir + "contacts_" + str(next_time_step) + ".dou", 'w') as fp:
                    fp.write(''.join((one_file_reader(
                        simulation_directory + '/contacts/contacts_' + str(next_time_step) + '.dou'))))

        else: print("Contacts not read as they are not to be plotted")

        if (one_file_reader(
                 simulation_directory + '/mirror_particles/mirror_particles_' + str(time_step) + '.dou') == OSError):
            print("file not existing")
        else:
             # =========WRITING MIRROR PARTICLES========================================================================
            print("Writing mirror particle file")
            with open(mirror_particles_dir + 'mirror_particles_' + str(time_step) + '.dou', 'w') as fp:
                fp.write(''.join((one_file_reader(
                    simulation_directory + '/mirror_particles/mirror_particles_' + str(time_step) + '.dou'))))

            # Write previous particle_time
            if not os.path.exists(mirror_particles_dir + "mirror_particles_" + str(previous_time_step) + ".dou"):
                with open(mirror_particles_dir + "mirror_particles_" + str(previous_time_step) + ".dou", 'w') as fp:
                    fp.write(''.join((one_file_reader(
                        simulation_directory + '/mirror_particles/mirror_particles_' + str(
                            previous_time_step) + '.dou'))))

            # Write next particle_time
            if not os.path.exists(mirror_particles_dir + "mirror_particles_" + str(next_time_step) + ".dou"):
                with open(mirror_particles_dir + "mirror_particles_" + str(next_time_step) + ".dou", 'w') as fp:
                    fp.write(''.join((one_file_reader(
                        simulation_directory + '/mirror_particles/mirror_particles_' + str(
                            next_time_step) + '.dou'))))

            # =========WRITING PERIODIC BC==============================================================================
            if not os.path.exists(local_directory + "periodic_bc.dou"):
                print("Writing periodic bc file")
                with open(local_directory + "periodic_bc.dou", 'w') as fp:
                    fp.write(''.join((one_file_reader(simulation_directory + '/periodic_bc.dou'))))
            else:
                print('periodic_BC already local')



    # =========COMPLEMENTING PARTICLE CONTACTS==========================================================================
    elif plot_contacts and not os.path.exists(contacts_dir + "contacts_" + str(time_step) + ".dou"):
        print("Complementing temporary files with contacts files")
        with open(contacts_dir + "contacts_" + str(time_step) + ".dou", 'w') as fp:
            fp.write(''.join((one_file_reader(simulation_directory + '/contacts/contacts_' + str(time_step) + '.dou'))))
    else: print('Local directory exists: '+local_directory)
    # ==================================================================================================================

    fig = mlab.figure(size=(1920, 1200), bgcolor=(1., 1., 1.), fgcolor=(0, 0., 0.))
    visual.set_viewer(fig)

    if plot_contacts:
        snapshot = Snapshot(os.path.expanduser(local_directory),BatteryContactPlotter)
        snapshot.contact_plotter.color = colors.white
        snapshot.contact_plotter.binder_radius = 0.02
    else: snapshot = Snapshot(os.path.expanduser(local_directory))

    snapshot.particle_opacity = 1
    snapshot.mirror_particles_plotter.opacity = 1
    snapshot.mirror_particles = mirror_particles
    snapshot.plot_periodic_bc=periodic_BCs
    snapshot.plot_velocity_arrow=velocity_arrow
    snapshot.plot_force_arrow= force_arrow
    # snapshot.surfaces_colors[1] = colors.red
    # snapshot.surfaces_opacities[1] = 0
    # snapshot.surfaces_opacities[0] = 0
    #Initiate the bounding box
    bbox = BoundingBox()
    #Set limits in the x,y and z direction
    bbox.x_max = lambda t: 2.5
    bbox.x_min = lambda t: -2.5
    bbox.y_max = lambda t: 2.5
    bbox.y_min = lambda t: -2.5
    # bbox.z_max = lambda t: 2
    #Apply bounding box to surfaces
    snapshot.surface_bounding_boxes[0] = bbox
    snapshot.surface_bounding_boxes[1] = bbox
    snapshot.surface_bounding_boxes[2] = bbox
    snapshot.surface_bounding_boxes[3] = bbox
    snapshot.surface_bounding_boxes[4] = bbox
    snapshot.surface_bounding_boxes[5] = bbox
    snapshot.plot(time_step)

    # ==SHOWING PLOT====================================================================================================
    print('showing plot')
    mlab.show()
