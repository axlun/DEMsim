//
// Created by Axel on 2022-05-18.
//

#include "../simulations.h"

#include <vector>

#include "../../engine/engine.h"
#include "../../contact_models/viscoelastic.h"
#include "../../contact_models/Binder_behavour_investigation/elastic_plastic_binder_rigid_perfect_plastic_particle.h"
#include "../../materials/electrode_material.h"


void DEM::restart_test(const std::string &settings_file_name)
{
    using namespace DEM;
    using ForceModel = elastic_plastic_binder_rigid_perfect_plastic_particle;
    using ParticleType = SphericalParticle<ForceModel>;
    using EngineType = Engine<ForceModel,ParticleType>;
    using namespace std::chrono_literals;

    SimulationParameters parameters(settings_file_name);
    auto restart_file_name =parameters.get_parameter<std::string>("restart_file_name");
    auto output_directory = parameters.get_parameter<std::string>("output_dir");
    auto simulator = EngineType(restart_file_name);
    auto mat = simulator.create_material<ElectrodeMaterial>(4800);
    auto Calendering_output = simulator.get_output("output_0");
    simulator.remove_output(Calendering_output);
    auto restart_output = simulator.create_output(output_directory+"restart_file", 0.05s);
    restart_output->print_particles = true;
    restart_output->print_surface_positions = true;
    restart_output->print_kinetic_energy = true;
    restart_output->print_contacts = true;
    restart_output->print_surface_forces = true;
    restart_output->print_fabric_force_tensor =true;
    restart_output->print_periodic_bc = true;
    restart_output->print_mirror_particles= true;

    auto top_surface = simulator.get_surface<EngineType::PointSurfacePointer>("top_plate");
    auto deformable_surface = simulator.get_surface<EngineType::DeformablePointSurfacePointer>("bottom_plate");
    std::cout << "****************Time to start the simulation**************** \n";

    EngineType::RunForTime run_time(simulator, 1.0s);
    simulator.run(run_time);

}