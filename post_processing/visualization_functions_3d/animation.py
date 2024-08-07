from collections import defaultdict
import glob
import os
import re
import time
import pathlib

import numpy as np
from mayavi import mlab

from visualization_functions_3d.plotting_functions import SpheresPlotter, BoundingBox, SurfacesPlotter, \
    ForceArrowPlotter, ForceArrowPlotter2, VelocityArrowPlotter
from visualization_functions_3d.periodic_bc import PeriodicBC
from visualization_functions_3d import colors


class Animation:
    def __init__(self, directory):
        self.directory = pathlib.Path(directory)
        self.delay = 0.0
        self.nth_timeframe = 1
        self.start_time = 0
        self.end_time = None
        self.max_decimal_points = 0
        self.save_frames = False
        self.frame_times = None
        self.save_directory = ''
        self.image_file_prefix = 'frame'
        self.image_file_extension = 'png'
        self.figure_directory = None
        self.surfaces_colors = defaultdict(lambda: (0., 0., 1.))
        self.surfaces_opacities = defaultdict(lambda: 0.5)
        self.plot_order = None
        self.bounding_boxes = defaultdict(BoundingBox)
        self.visible_functions = defaultdict(lambda: lambda t: True)
        self.dpi = 500, 800
        self.zoom_settings = None
        self.plot_periodic_bc = False
        self.periodic_bc_plotter = None
        self.force_arrow_plotter = None
        self.velocity_arrow_plotter = None
        self.mirror_particles = False
        self.plot_force_arrow = False
        self.plot_velocity_arrow = False
        self.n_vel = -1
        self.n_force = -1
        self.spheres_plotter = SpheresPlotter()
        self.mirror_particles_plotter = SpheresPlotter(color=colors.silver)

        self.initialized = False
        self.surfaces_plotter = None
        self.view_surfaces = True

    def run(self):
        a = self._animation()
        for _ in a:
            pass

    def initialize(self):
        if self.view_surfaces:
            self.surfaces_plotter = SurfacesPlotter(self.directory / 'surface_positions.dou', self.surfaces_colors,
                                                    self.surfaces_opacities, self.plot_order, self.bounding_boxes,
                                                    self.visible_functions)
            self.surfaces_plotter.surfaces_opacities = self.surfaces_opacities
            self.surfaces_plotter.surfaces_colors = self.surfaces_colors
            self.surfaces_plotter.plot_order = self.plot_order
            self.surfaces_plotter.bounding_boxes = self.bounding_boxes
            self.surfaces_plotter.visible_times = self.visible_functions

        if self.plot_periodic_bc:
            self.periodic_bc_plotter = PeriodicBC(self.directory / 'periodic_bc.dou')

        if self.plot_force_arrow:
            # self.force_arrow_plotter = ForceArrowPlotter(colors.green, 1., n_force=self.n_force)
            self.force_arrow_plotter = ForceArrowPlotter2(colors.green, 1.,)

        if self.plot_velocity_arrow:
            self.velocity_arrow_plotter = VelocityArrowPlotter(colors.red, 1., n_velocity=self.n_vel)

        particle_files = glob.glob(str(self.directory) + '/particles/particles_*.dou')
        particle_files = [os.path.basename(particle_file) for particle_file in particle_files]
        print(self.nth_timeframe)
        particle_files = particle_files[::self.nth_timeframe]

        # Getting frametimes
        self.frame_times = []
        for p_file in particle_files:
            #           self.frame_times.append(re.findall(r"[-+]?\d*\.\d+|\d+", p_file)[0])
            time_frame_split = re.split(r"particles_|.dou", p_file)[1]
            self.frame_times.append(time_frame_split)
            if time_frame_split[::-1].find('.') > self.max_decimal_points:
                self.max_decimal_points = time_frame_split[::-1].find('.')

        self.frame_times = np.array(sorted(self.frame_times, key=lambda x: float(x)), dtype=str)
        frame_times_np = np.array([float(t) for t in self.frame_times])
        self.frame_times = self.frame_times[frame_times_np >= self.start_time]
        frame_times_np = frame_times_np[frame_times_np >= self.start_time]
        if self.end_time:
            self.frame_times = self.frame_times[frame_times_np < self.end_time]
        print(self.frame_times)

        if self.save_frames and not os.path.isdir(self.save_directory):
            os.makedirs(self.save_directory)
        self.initialized = True

    def _animation(self):
        if not self.initialized:
            self.initialize()
        n = len(self.frame_times)

        file_name_str = []
        for i, t in enumerate(self.frame_times):
            string = str(t)
            if string[::-1].find('.') <= 0:
                string += ".0"
            if string[::-1].find('.') < self.max_decimal_points:
                string += "0"*(self.max_decimal_points-int(string[::-1].find('.')))
            file_name_str.append(string)

        surface_position_data = np.array([])
        if self.plot_velocity_arrow:
            surface_position_data = np.genfromtxt(self.directory / 'surface_positions.dou', delimiter=',')
        for i, t in enumerate(self.frame_times):
            print('Frame number: '+ str(i))
            print('Frame time: ' + str(t))
            particle_data = np.genfromtxt(self.directory / str('particles/particles_' + str(t) + '.dou'), delimiter=',')
            self.spheres_plotter.plot(particle_data)
            if self.mirror_particles:
                mirror_particle_data = np.genfromtxt(self.directory / str('mirror_particles/mirror_particles_'
                                                     + t + '.dou'), delimiter=',')
                self.mirror_particles_plotter.plot(mirror_particle_data)
            if self.surfaces_plotter:
                self.surfaces_plotter.plot(float(t))
            if self.plot_periodic_bc:
                self.periodic_bc_plotter.plot(float(t))
            if self.plot_force_arrow:
                self.force_arrow_plotter.plot(particle_data)
            if self.plot_velocity_arrow:
                self.velocity_arrow_plotter.plot(self.directory, float(t), particle_data, surface_position_data)

            f = mlab.gcf()
            f.scene.render()
            # mlab.colorbar()
            # mlab.axes()
            if self.save_frames:
                if self.figure_directory is None:
                    self.figure_directory = self.directory
                if not os.path.isdir(self.figure_directory):
                    os.makedirs(self.figure_directory)

                #each file name is for the frame time in the simulation
                # name = '/' + self.image_file_prefix + str(self.frame_times[i]).replace(".","-") + '.' + self.image_file_extension
                name = '/' + self.image_file_prefix + file_name_str[i] + '.' + self.image_file_extension
                # each file has name frame_00x, _0xx, xxx etc
                # name = '/' + self.image_file_prefix + '0'*(len(str(n))-len(str(i))) + str(i) + '.' \
                #        + self.image_file_extension
                mlab.savefig(filename=self.save_directory + name)

            time.sleep(self.delay)
            if self.plot_velocity_arrow:
                self.velocity_arrow_plotter.clear_arrows()
            # if self.plot_force_arrow:
            #     self.force_arrow_plotter.clear_arrows()


            yield