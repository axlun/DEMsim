from mayavi import mlab

from visualization_functions_3d.animation import Animation

def main():
    print('start')
#    simulation_directory = 'c:/Users/Axel/Documents/DEM/results/electrode_calendaring/restart_test/restart_file'
#    simulation_directory = 'c:/Users/Axel/Documents/DEM/results/electrode_calendaring/SN0_250particlesWithBinder_new_packingmethod_with_periodic_BC'
#    simulation_directory = 'c:/Users/Axel/Documents/DEM/results/electrode_calendaring/SN0_1000particlesWithBinder_new_packingmethod_2'
#    simulation_directory = 'c:/Users/Axel/Desktop/SN00_400p_plastic_binder'
    simulation_directory = 'c:/Users/Axel/Documents/DEM/results/electrode_mechanical_loading/SN00_test'
    mlab.figure(size=(1920, 1200), bgcolor=(1., 1., 1.))
    animation = Animation(simulation_directory)
    animation.save_directory = simulation_directory+'/animation/imgs/'
    animation.save_frames = True
    animation.delay = 1e-6
    animation.plot_periodic_bc = True
    animation.mirror_particles = False
    animation.view_surfaces = True
    animation.start_time = 14.4
#    animation.end_time = animation.start_time +.3
    animation.run()
    mlab.show()

if __name__ == '__main__':
    main()