from Bertil_functions import *



def simulation_runner(program, sim_file):



    input_command = "cd DEMsim/build\n./rebuild.sh\nqsub -v prog="+program+",simulation="+sim_file+" run_dem_simulation.sh\n"
    commad_input(input_command)

    return
if __name__ == '__main__':

    simulation_name = "electrode_calendering"
    simulation_file = "/DEMsim/simulations/force_model_impact_on_electrode/Electrode_calendering_bertil"#dont add .sim
    #simulation_name = "electrode_mechanical_loading"
    #simulation_file = "/DEMsim/simulations/force_model_impact_on_electrode/Electrode_mechanical_loadingYes_bertil"  # dont add .sim

    simulation_file = "/scratch/users/axlun/"+simulation_file
    simulation_runner(simulation_name,simulation_file)


    # command0 = "cd /scratch/users/axlun/DEMsim/"
    # command = "ls"
    # input_command = command0 +" \n" + command
    # input_command = "cd DEMsim/build\n./rebuild.sh\nls"
    # commad_input(input_command)