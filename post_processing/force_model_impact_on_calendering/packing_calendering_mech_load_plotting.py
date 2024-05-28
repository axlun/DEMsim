from force_model_impact_on_calendering.Bertil_calendering_pressure import calendering_pressure_plotting_func
from force_model_impact_on_calendering.Bertil_packing_process import packing_process_plotting_func
from force_model_impact_on_calendering.Bertil_mechanical_properties import mechanical_properties_plotting_func
import matplotlib.pyplot as plt






if __name__=='__main__':
    # simulation_directory = ('/scratch/users/axlun/DEMsim/results/article_2/final_runs_2/SN_101/1/','hertz')
    # simulation_directory = ('/scratch/users/axlun/DEMsim/results/article_2/final_runs_2/SN_101/2/','hertz')
    # simulation_directory = ('/scratch/users/axlun/DEMsim/results/article_2/final_runs_2/SN_101/3/','hertz')
    # simulation_directory = ('/scratch/users/axlun/DEMsim/results/article_2/final_runs_2/SN_101/4/','hertz')

    # simulation_directory = ('/scratch/users/axlun/DEMsim/results/article_2/final_runs_2/SN_111/1/','hertz')
    # simulation_directory = ('/scratch/users/axlun/DEMsim/results/article_2/final_runs_2/SN_111/2/','hertz')
    # simulation_directory = ('/scratch/users/axlun/DEMsim/results/article_2/final_runs_2/SN_111/3/','hertz')
    # simulation_directory = ('/scratch/users/axlun/DEMsim/results/article_2/final_runs_2/SN_111/4/','hertz')

    # simulation_directory = ('/scratch/users/axlun/DEMsim/results/article_2/final_runs_2/SN_301/1/',
    #                         'el_pl_binder_el_pl_particle')
    # simulation_directory = ('/scratch/users/axlun/DEMsim/results/article_2/final_runs_2/SN_301/2/',
    #                         'el_pl_binder_el_pl_particle')
    # simulation_directory = ('/scratch/users/axlun/DEMsim/results/article_2/final_runs_2/SN_301/3/',
    #                         'el_pl_binder_el_pl_particle')
    # simulation_directory = ('/scratch/users/axlun/DEMsim/results/article_2/final_runs_2/SN_301/4/',
    #                         'el_pl_binder_el_pl_particle')
    # simulation_directory = ('/scratch/users/axlun/DEMsim/results/article_2/final_runs_2/SN_301/5/',
    #                         'el_pl_binder_el_pl_particle')
    # simulation_directory = ('/scratch/users/axlun/DEMsim/results/article_2/final_runs_2/SN_301/6/',
    #                         'el_pl_binder_el_pl_particle')

    # simulation_directory = ('/scratch/users/axlun/DEMsim/results/article_2/final_runs_2/SN_201/1/',
    #                         'el_pl_binder_el_pl_particle')
    # simulation_directory = ('/scratch/users/axlun/DEMsim/results/article_2/final_runs_2/SN_201/2/',
    #                         'el_pl_binder_el_pl_particle')
    # simulation_directory = ('/scratch/users/axlun/DEMsim/results/article_2/final_runs_2/SN_201/3/',
    #                         'el_pl_binder_el_pl_particle')
    # simulation_directory = ('/scratch/users/axlun/DEMsim/results/article_2/final_runs_2/SN_201/4/',
    #                         'el_pl_binder_el_pl_particle')

    # simulation_directory = ('/scratch/users/axlun/DEMsim/results/article_2/final_runs_2/SN_200/1/',
    #                         'el_pl_binder_el_pl_particle')
    # simulation_directory = ('/scratch/users/axlun/DEMsim/results/article_2/final_runs_2/SN_200/2/',
    #                         'el_pl_binder_el_pl_particle')
    # simulation_directory = ('/scratch/users/axlun/DEMsim/results/article_2/final_runs_2/SN_200/3/',
    #                         'el_pl_binder_el_pl_particle')
    # simulation_directory = ('/scratch/users/axlun/DEMsim/results/article_2/final_runs_2/SN_200/4/',
    #                         'el_pl_binder_el_pl_particle')

    # simulation_directory = ('/scratch/users/axlun/DEMsim/results/article_2/final_runs_2/SN_202/1/',
    #                         'el_pl_binder_el_pl_particle')
    # simulation_directory = ('/scratch/users/axlun/DEMsim/results/article_2/final_runs_2/SN_202/2/',
    #                         'el_pl_binder_el_pl_particle')
    # simulation_directory = ('/scratch/users/axlun/DEMsim/results/article_2/final_runs_2/SN_202/3/',
    #                         'el_pl_binder_el_pl_particle')
    # simulation_directory = ('/scratch/users/axlun/DEMsim/results/article_2/final_runs_2/SN_202/4/',
    #                         'el_pl_binder_el_pl_particle')

    # simulation_directory = ('/scratch/users/axlun/DEMsim/results/article_2/final_runs_2/SN_101P/1/',
    #                         'hertz')
    # simulation_directory = ('/scratch/users/axlun/DEMsim/results/article_2/final_runs_2/SN_101P/2/',
    #                         'hertz')
    # simulation_directory = ('/scratch/users/axlun/DEMsim/results/article_2/final_runs_2/SN_101P/3/',
    #                         'hertz')
    # simulation_directory = ('/scratch/users/axlun/DEMsim/results/article_2/final_runs_2/SN_101P/4/',
    #                         'hertz')



    # simulation_directory = ('/scratch/users/axlun/DEMsim/results/article_2/final_runs_2/SN_2011/1/',
    #                         'el_pl_binder_el_pl_particle')
    # simulation_directory = ('/scratch/users/axlun/DEMsim/results/article_2/final_runs_2/SN_2011/2/',
    #                         'el_pl_binder_el_pl_particle')
    # simulation_directory = ('/scratch/users/axlun/DEMsim/results/article_2/final_runs_2/SN_2011/3/',
    #                         'el_pl_binder_el_pl_particle')
    # simulation_directory = ('/scratch/users/axlun/DEMsim/results/article_2/final_runs_2/SN_2011/4/',
    #                         'el_pl_binder_el_pl_particle')


    # simulation_directory = ('/scratch/users/axlun/DEMsim/results/article_2/final_runs_2/SN_2111/1/',
    #                         'el_pl_binder_el_pl_particle')
    # simulation_directory = ('/scratch/users/axlun/DEMsim/results/article_2/final_runs_2/SN_2111/2/',
    #                         'el_pl_binder_el_pl_particle')
    # simulation_directory = ('/scratch/users/axlun/DEMsim/results/article_2/final_runs_2/SN_2111/3/',
    #                         'el_pl_binder_el_pl_particle')
    simulation_directory = ('/scratch/users/axlun/DEMsim/results/article_2/final_runs_2/SN_2111/4/',
                            'el_pl_binder_el_pl_particle')


    packing_process_plotting_func(simulation_directory[0] + 'electrode_natural_packing_'+simulation_directory[1]+'/')
    calendering_pressure_plotting_func(simulation_directory[0] + 'electrode_calendering_'+simulation_directory[1]+'/')
    # mechanical_properties_plotting_func(simulation_directory[0] + 'electrode_mechanical_loading_' +
    #                                     simulation_directory[1])
    plt.show()