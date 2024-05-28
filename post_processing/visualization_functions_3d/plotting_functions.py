import statistics
import sys

import matplotlib.pyplot as plt

from collections import namedtuple
from collections import OrderedDict
from math import pi

import numpy as np

from mayavi import mlab
from tvtk.tools import visual

from visualization_functions_3d import colors
import time


# todo create funciton to move arrows instead of creating new each timestep
def arrow_creator(p1, v_p1, color=colors.red, radius=1e-1,  cone_radius=2e-1,scale_factor=1.0):
    # start_time = time.time()
    p1 = np.array(p1)
    # print('p1:' + str(p1))
    v_p1 = np.array(v_p1)
    # print('V_p1:' + str(v_p1))
    # create_time = time.time()
    ar1 = visual.arrow(x=p1[0], y=p1[1], z=p1[2])
    # print('Create time = '+str(time.time()-create_time)+'s \n')
    arrow_length = scale_factor*np.linalg.norm(v_p1)
    # print('Arrow length:' + str(arrow_length))
    ar1.length_cone = cone_radius
    ar1.radius_shaft = radius
    ar1.radius_cone = cone_radius
    if arrow_length != 0:
        ar1.actor.scale = [arrow_length, arrow_length, arrow_length]
        ar1.pos = ar1.pos/arrow_length
    ar1.axis = v_p1
    ar1.color = color
    # print('Execute time = '+str(time.time()-start_time)+'s \n')
    return ar1


class BoundingBox:
    def __init__(self):
        self.x_min = lambda t: -1e99
        self.x_max = lambda t: 1e99
        self.y_min = lambda t: -1e99
        self.y_max = lambda t: 1e99
        self.z_min = lambda t: -1e99
        self.z_max = lambda t: 1e99

    def values(self, time):
        return [self.x_min(time), self.x_max(time),
                self.y_min(time), self.y_max(time),
                self.z_min(time), self.z_max(time)]


# Todo, fix unused parameter opacity
class SpheresPlotter:
    def __init__(self, opacity=1., color=colors.copper):
        self.ms = None
        self.color = color
        self.opacity = opacity

    def plot(self, data, function=None, vmin=None, vmax=None):
        if len(data) > 0:
            if len(data.shape) == 1:
                data = np.expand_dims(data, 0)
            x = data[:, 1]
            y = data[:, 2]
            z = data[:, 3]
            r = data[:, 7]
            if vmin is None:
                vmin = np.min(function)
            if vmax is None:
                vmax = np.max(function)
            if function is not None:
                pts = mlab.points3d(x, y, z, resolution=32, opacity=self.opacity, transparent=self.opacity != 1.,
                                    scale_factor=1., vmax=vmax, vmin=vmin)
                self.ms = pts.mlab_source
                pts.glyph.scale_mode = 'scale_by_vector'
                self.ms.dataset.point_data.vectors = np.tile(2/np.sqrt(3)*r, (3, 1)).transpose()
                self.ms.dataset.point_data.scalars = function

            else:
                if self.ms is None: # If no point -> Create points
                    self.ms = mlab.points3d(x, y, z, 2*r,
                                            color=self.color,
                                            resolution=32,
                                            scale_factor=1.,
                                            scale_mode='scalar',
                                            opacity=self.opacity,
                                            reset_zoom=False).mlab_source
                elif self.ms.points.shape[0] == data.shape[0]: # If points -> update x,y,z
                    self.ms.set(x=x, y=y, z=z, scalars=2*r)
                else:
                    self.ms.reset(x=x, y=y, z=z, scalars=2*r)


def fulfill_bounding_box(bounding_box, x, y, z, time):
    values = bounding_box.values(time)
    x[x < values[0]] = values[0]
    x[x > values[1]] = values[1]
    y[y < values[2]] = values[2]
    y[y > values[3]] = values[3]
    z[z < values[4]] = values[4]
    z[z > values[5]] = values[5]


class PointSurfacePlotter:
    def __init__(self, bounding_box=None):
        self.ms = None
        self.bounding_box = bounding_box

    def plot(self, data, color, opacity, time=0.):
        n = int(data.shape[0]/3)
        x = data[0:3*n-2:3]
        y = data[1:3*n-1:3]
        z = data[2:3*n:3]
        if self.bounding_box:
            fulfill_bounding_box(self.bounding_box, x, y, z, time)

        # This orders the points so that a rectangle is plotted
        # pts = mlab.points3d(x, y, z)
        # mesh = mlab.pipeline.delaunay2d(pts)
        #
        # pts.remove()

        # =SURFACE PLOTTING WITH TRIANGULAR MESH========================================================================
        # Creates two two trianges for a surface with four vertices,
        triangles = [[0,1,2],[2,3,0]] # verticies index for each triangle

        if self.ms is None:
            self.ms = mlab.triangular_mesh(x,y,z,triangles,color=color,opacity=opacity,reset_zoom=False)

            # mlab.triangular_mesh(x, y, z, triangles, color=(0,0,0),representation='wireframe' , reset_zoom=False,line_width=5)
            # ^ use to also plot wireframe for triangular mesh
            # self.ms = mlab.pipeline.surface(mesh, color=color, opacity=opacity, transparent=True,
            #                                 reset_zoom=False).mlab_source
        else:
            # Updating the pipeline with the new set of points
            # There is probably a cuter way to do this
            # self.ms.points = mesh.mlab_source.points


            self.ms.mlab_source.x = x
            self.ms.mlab_source.y = y
            self.ms.mlab_source.z = z
            # Update the x,y,z coordinates for the verticies, connectivity is not changed

class CylinderPlotter:
    def __init__(self, bounding_box=None):
        self.ms = None
        self.upper_plate = None
        self.lower_plate = None
        self.bounding_box = bounding_box
        self.length_extension = 0.
        self.closed = False

    def plot(self, data, color, opacity, time=0.):
        r = data[0]
        # axis = data[1:4]   # Todo use axis parameter
        point = data[4:7]
        length = data[7] + self.length_extension

        q, z = np.meshgrid(np.linspace(0, 2*pi, 100), np.linspace(point[2], point[2] + length, 100))

        x = r*np.cos(q) + point[0]
        y = r*np.sin(q) + point[1]
        if self.bounding_box:
            fulfill_bounding_box(self.bounding_box, x, y, z, time)
        if self.ms is None:
            self.ms = mlab.mesh(x, y, z, color=color, opacity=opacity, transparent=True, reset_zoom=False).mlab_source
        else:
            self.ms.set(x=x, y=y, z=z)

        if self.closed:
            rad, q = np.meshgrid(np.linspace(0, r, 100), np.linspace(0, 2*pi, 100))
            x = rad*np.cos(q) + point[0]
            y = rad*np.sin(q) + point[1]

            if self.upper_plate is None:
                self.upper_plate = mlab.mesh(x, y, 0*x + point[2], color=color, opacity=opacity, transparent=True,
                                             reset_zoom=False).mlab_source
            else:
                self.upper_plate.set(x=x, y=y, z=z)

            if self.lower_plate is None:
                self.lower_plate = mlab.mesh(x, y, 0*x, color=color, opacity=opacity, transparent=True,
                                             reset_zoom=False).mlab_source
            else:
                self.lower_plate.set(x=x, y=y, z=z)


PlotObject_ = namedtuple('PlotObject', ['start_idx', 'end_idx'])


class SurfacesPlotter:
    def __init__(self, surface_file_name, surfaces_colors=None, surfaces_opacities=None, plot_order=None,
                 bounding_boxes=None, visible_times=None):

        self.plotters = {}
        self.plotter_data = {}
        self.data = OrderedDict()
        self.counter = 0
        self.set_data_file(surface_file_name)

        self.surfaces_colors = {}
        self.surfaces_opacities = {}
        self.plot_order = plot_order
        self.bounding_boxes = {}
        self.visible_times = {}

        if surfaces_colors is None:
            surfaces_colors = {}
        if surfaces_opacities is None:
            surfaces_opacities = {}
        if bounding_boxes is None:
            bounding_boxes = {}
        if visible_times is None:
            visible_times = {}

        for surface_id in self.plotters:
            self.surfaces_colors[surface_id] = surfaces_colors.get(surface_id, (0., 0., 1.))
            self.surfaces_opacities[surface_id] = surfaces_opacities.get(surface_id, 0.5)
            self.bounding_boxes[surface_id] = bounding_boxes.get(surface_id, BoundingBox())
            self.visible_times[surface_id] = visible_times.get(surface_id, lambda t: True)

    def set_data_file(self, surface_file_name):
        with open(surface_file_name) as data_file:
            data_lines = data_file.readlines()

        # Inspect the first line
        line = data_lines[0]
        words = line.split(", ")
        # print(words)
        id_idx = [i for i in range(len(words)) if words[i].upper().startswith('ID')]
        for idx in id_idx:
            surface_id = int(words[idx][3:])
            surface_type = words[idx+1][5:]

            if surface_type == 'Cylinder':
                self.plotter_data[surface_id] = PlotObject_(idx+2, idx+10)
                self.plotters[surface_id] = CylinderPlotter()
            elif surface_type in ["PointSurface", "DeformablePointSurface"]:
                num_points = int(words[idx+2])
                self.plotter_data[surface_id] = PlotObject_(idx+3, idx+3+num_points*3)
                self.plotters[surface_id] = PointSurfacePlotter()

        data = np.genfromtxt(surface_file_name, delimiter=',')
        for i in range(data.shape[0]):
            if len(data.shape) == 1:
                self.data[data[i]] = data[i]
            else:
                self.data[data[i, -1]] = data[i, :-1]

    def plot(self, t=None):
        if t is None:
            data_line = self.data[self.data.keys()[self.counter]]
            self.counter += 1
        else:
            try:
                data_line = self.data[t]
                self.counter = 0
            except LookupError:
                print("No surface data at time ", t)
                sys.exit(1)

        plot_order = self.plot_order
        if self.plot_order is None:
            plot_order = self.plotters.keys()

        for surface_id in plot_order:
            if self.visible_times[surface_id](t):
                plotter = self.plotters[surface_id]
                plotter.bounding_box = self.bounding_boxes[surface_id]
                plotter_data = self.plotter_data[surface_id]
                data = data_line[plotter_data.start_idx:plotter_data.end_idx]
                plotter.plot(data, self.surfaces_colors[surface_id], self.surfaces_opacities[surface_id], t)


class VelocityArrowPlotter:
    def __init__(self, color=colors.green, opacity=1.,n_velocity=-1):
        self.ms = None
        self.color = color
        self.opacity = opacity
        self.n_velocity = n_velocity
        self.velocity_arrow_array = []

    def plot(self,directory,time,particle_data,surface_position_data):
        # Calculate velocity vector from t-dt and t+dt
        #Find time index in surface position data and extract the time before and after
        time_step_index = np.where(surface_position_data[:,-1] == time)[0][0]

        side_length = abs(2 * surface_position_data[0,6])
        prev_time_step_index = time_step_index - 1
        next_time_step_index = time_step_index + 1
        if prev_time_step_index < 0:
            prev_time_step_index = time_step_index

        if next_time_step_index > np.size(surface_position_data[:,-1])-1:
            next_time_step_index = [np.size(surface_position_data[:,-1])-1]
        #calculateing time increment
        next_time_step = surface_position_data[next_time_step_index,-1]
        previous_time_step = surface_position_data[prev_time_step_index,-1]
        dt = next_time_step - previous_time_step

        if next_time_step.is_integer():
            next_time_step = int(next_time_step)
        if previous_time_step.is_integer():
            previous_time_step = int(previous_time_step)
        print('current timestep = ' + str(time))
        print('Next timestep = ' + str(next_time_step))
        print('Previous timestep = ' + str(previous_time_step))

        a = str('particles/particles_' + str(next_time_step) + '.dou')
        b = str('particles/particles_' + str(previous_time_step) + '.dou')

        next_particle_data = np.genfromtxt(directory / a , delimiter=',')
        previous_particle_data = np.genfromtxt(directory / b, delimiter=',')

        x = particle_data[:, 1]
        y = particle_data[:, 2]
        z = particle_data[:, 3]
        V_x = np.zeros(x.size)
        V_y = np.zeros(x.size)
        V_z = np.zeros(x.size)


        scale_factor = 1
        for i in range(0,x.size):
# ====================X=================================================================================================
            V_x[i] = scale_factor*Periodic_bc_velocity_correction(previous_particle_data[i, 1],next_particle_data[i, 1],side_length,dt)
# ======================================================================================================================

#======================Y================================================================================================
            V_y[i] = scale_factor*Periodic_bc_velocity_correction(previous_particle_data[i, 2], next_particle_data[i, 2], side_length, dt)
#=======================================================================================================================

#======================Z================================================================================================
            V_z[i] = scale_factor*(next_particle_data[i, 3] - previous_particle_data[i, 3]) / dt
#=======================================================================================================================

        V_length = (V_x ** 2 + V_y ** 2 + V_z ** 2) ** 0.5
        V_max = max(V_length)
        print('\nV_max = '+str(V_max))
        V_counter = 0

        V_n_th = np.sort(V_length)[-self.n_velocity]
        for i in range(0,x.size):
            if V_length[i] >= V_max/2 and V_length[i] >= V_n_th:
            # if particle_data[i,10] == 0:
            # if particle_data[i,0] in (236, 417):
            #if particle_data[i,0] == 417:
                print("Particle ID: " + str(int(particle_data[i,0])) + ' with v=' + str(V_length[i]))
                V_counter += 1
                self.velocity_arrow_array.append(arrow_creator([x[i], y[i], z[i]],  [V_x[i], V_y[i], V_z[i]]))
        print('Number of velocity vectors plotted:'+str(V_counter)+'\n')
        # velocity_arrow()

    def clear_arrows(self):
        for i in self.velocity_arrow_array:
            visual.remove_actor(i)


class ForceArrowPlotter2:
    def __init__(self,color=colors.red, opacity=1):
        self.ms = None
        self.color = color
        self.opacity = opacity

    def plot(self, data, function=None, vmin=None, vmax=None):
        if len(data) > 0:
            if len(data.shape) == 1:
                data = np.expand_dims(data,0)
        x = data[:, 1]
        y = data[:, 2]
        z = data[:, 3]
        f_x = data[:, 10]
        f_y = data[:, 11]
        f_z = data[:, 12]
        # ====================
        # Rescale force here
        F_length = (f_x**2+f_y**2+f_z**2)**0.5
        max_F =np.amax(F_length)
        max_particle_rad = np.amax(data[:,7])
        scale_factor = 10
        arrow_scale_factor = scale_factor * max_particle_rad/max_F

        f_x = f_x*arrow_scale_factor
        f_y = f_y*arrow_scale_factor
        f_z = f_z*arrow_scale_factor

        # print("max_F: " + str(max_F))
        # print("max_particle_R: " + str(max_particle_rad))
        # print('arrow_scale_factor: ' + str(arrow_scale_factor))
        # print('Longest arrow: ' + str(arrow_scale_factor*max_F))
        # ====================
        if vmin is None:
            vmin = np.min(function)
        if vmax is None:
            vmax = np.min(function)
        if function is not None: # to use if there is a vector function u,v,w(x,y,z)
            print('Not implemented')
            # arrw = mlab.quiver3d(x,y,z,f_x,f_y,f_z, resolution=32,opacity=self.opacity,
            #                      transparency=self.opacity != 1., scale_factor=1., vmax=vmax,vmin=vmin, scale_mode='vector')
            # self.ms = arrw.mlab_source
            # arrw.glyph.scale_mode = 'scale_by_vectors'
            # self.ms.dataset.point_data.vectors = np.tile(2 / np.sqrt(3) * r, (3, 1)).transpose()    # Fix these lines
            # self.ms.dataset.point_data.scalars = function                                           #  i needed
        else:
            if self.ms is None: # If no arrows -> create arrows
                self.ms = mlab.quiver3d(x,y,z,f_x,f_y,f_z,
                                        mode='arrow',
                                        resolution=32,
                                        opacity=self.opacity, vmax=vmax, vmin=vmin, color=self.color,
                                        scale_factor=1.,
                                        scale_mode='vector',
                                        reset_zoom=False).mlab_source
            elif self.ms.points.shape[0] == data.shape[0]: # If arrow and they have same shape as next data point -> update X and U
                self.ms.set(x=x, y=y, z=z, u=f_x, v=f_y, w=f_z,scale_factor=1.)
            else:   # If arrows but if they have different shape, redraw
                self.ms.reset(x=x, y=y, z=z,u=f_x,v=f_y,w=f_z,scale_factor=1.)


class ForceArrowPlotter:
    def __init__(self, color=colors.green, opacity=1.,n_force = 1):
        self.ms = None
        self.color = color
        self.opacity = opacity
        self.n_force = n_force
        self.force_arrow_array = []

    def plot(self, particle_data):
        start_time = time.time()
        x = particle_data[:, 1]
        y = particle_data[:, 2]
        z = particle_data[:, 3]
        F_x = particle_data[:, 10]
        F_y = particle_data[:, 11]
        F_z = particle_data[:, 12]


        #Scale all vectors after the largest force and 5 times of the larges particle
        F_length = (F_x**2+F_y**2+F_z**2)**0.5
        median_F = statistics.median_high(F_length)
        max_F =np.amax(F_length)
        max_particle_rad = np.amax(particle_data[:,7])
        scale_factor = 10 * max_particle_rad/max_F

        print('\nF_median = '+ str(median_F))
        print('F_max = '+str(max_F))
        print('r_max = '+ str(max_particle_rad))
        print('Scale factor: '+str(scale_factor))

        F_x_scaled = F_x * scale_factor
        F_y_scaled = F_y * scale_factor
        F_z_scaled = F_z * scale_factor
        number_f_plotted = 0
        F_n_th = np.sort(F_length)[-self.n_force]


        # =PLOTTING OF ALL ARROWS=======================================================================================
        if len(self.force_arrow_array) == 0:
            for i in range(0,x.size):
                if F_length[i] >= max_F/10 and F_length[i] >= F_n_th:
                    number_f_plotted += 1

                    self.force_arrow_array.append(arrow_creator([x[i], y[i], z[i]], [F_x_scaled[i], F_y_scaled[i],
                                                                                     F_z_scaled[i]], color=colors.green))
            print('Number of force vectors plotted: ' + str(number_f_plotted)+'\n')
        # =MOVE ARROWS==============================================================================================
        else:
            print('Move arrows')
            # for i in range(0,x.size):
        # self.force_arrow_array.

    def clear_arrows(self):
        for i in self.force_arrow_array:
            visual.remove_actor(i)

def Periodic_bc_velocity_correction(prev_step,next_step,side_length,dt):
    # print(prev_step)
    if np.sign(prev_step) == 1 and np.sign(next_step) == -1 and abs(
            (next_step - prev_step)) > .9*side_length: # .9 * side_length means that particle travels very far in one time step
        V = ((next_step - prev_step) + side_length) / dt
        print('+ to - correction from ' + str((next_step - prev_step) / dt) + ' to ' + str(V))
        print('Previous step: ' + str(prev_step))
        print('Next step: ' + str(next_step))
    elif np.sign(prev_step) == -1 and np.sign(next_step) == 1 and abs(
            (next_step - prev_step)) > .9*side_length:
        V = ((next_step - prev_step) - side_length) / dt
        print('- to + correction from ' + str((next_step - prev_step) / dt) + ' to ' + str(V))
        print('Previous step: ' + str(prev_step))
        print('Next step: ' + str(next_step))

    else:
        V = (next_step - prev_step) / dt

    return V


if __name__ == '__main__':
    simulation_directory = '../results/cyclic_triaxial/test/'
    surfaces_plotter = SurfacesPlotter(simulation_directory + 'surface_positions.dat')
    bbox = BoundingBox()
    bbox.z_min = -0.01
    bbox.z_max = 0.05
    surfaces_plotter.plot()

    mlab.show()
