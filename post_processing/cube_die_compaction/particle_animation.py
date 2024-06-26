from mayavi import mlab

from visualization_functions_3d.animation import Animation


def main():
    simulation_directory = 'c:/Users/Axel/Documents/DEM/results/cube_die_compaction/3'
    mlab.figure(size=(1920, 1200), bgcolor=(1., 1., 1.))
    animation = Animation(simulation_directory)
    animation.save_directory = simulation_directory+'/animation/imgs/'
    animation.save_frames = True
    animation.delay = 0.0001
    # animation.plot_periodic_bc = True
    # animation.mirror_particles = False
    animation.start_time = 0
    animation.run()
    mlab.show()


if __name__ == '__main__':
    main()