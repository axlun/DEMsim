//
// Created by Axel on 2023-10-31.
//

#include "../../simulations.h"

#include <vector>

#include "../../../engine/engine.h"
#include "../../../contact_models/Positive_electrode/elastic_plastic_binder_elastic_plastic_particle.h"
#include "../../../materials/electrode_material.h"


void DEM::restart_file_tester(const std::string &settings_file_name)
{
    using namespace DEM;
    using ForceModel = elastic_plastic_binder_elastic_plastic_particle;
    using ParticleType = SphericalParticle<ForceModel>;
    using EngineType = Engine<ForceModel,ParticleType>;
    using namespace std::chrono_literals;

    SimulationParameters parameters(settings_file_name);

    auto restart_file_name =parameters.get_parameter<std::string>("restart_file_name");
    auto output_directory = parameters.get_parameter<std::string>("output_dir");
    std::cout << "restart_file_name = " << restart_file_name << std::endl;
    std::cout << "output_directory = " << output_directory << std::endl;
    auto simulator = EngineType(restart_file_name);

//    EngineType::RunForTime run_for_time(simulator, 1E-3s);
//    simulator.run(run_for_time);
    /*
    std::cout << "Write result file" << std::endl;
    std::ofstream results_file;
    results_file.open(output_directory + "/test_results_file.dou");
    results_file << "Test print in result file." << std::endl;
    results_file.close();
    */

    std::cout<<"Writing restart file " << std::endl;
    simulator.write_restart_file(output_directory + "/testing_restart_file.res");
    }