//
// Created by Axel on 2023-06-16.
//

#include "../../simulations.h"

#include <vector>

#include "../../../engine/engine.h"
#include "../../../contact_models/Positive_electrode/elastic_plastic_binder_elastic_plastic_particle.h"
#include "../../../materials/electrode_material.h"


void DEM::electrode_relaxation_el_pl_binder_el_pl_particle(const std::string &settings_file_name)
{
    using namespace DEM;
    using ForceModel = elastic_plastic_binder_elastic_plastic_particle;
    using ParticleType = SphericalParticle<ForceModel>;
    using EngineType = Engine<ForceModel,ParticleType>;
    using namespace std::chrono_literals;

    SimulationParameters parameters(settings_file_name);
    auto restart_file_name =parameters.get_parameter<std::string>("restart_file_name");
    auto output_directory = parameters.get_parameter<std::string>("output_dir");
    auto loading_direction = parameters.get_parameter<int>("loading_direction"); // compression->-1, tension->1
    auto strain_load = parameters.get_parameter<double>("strain");
    std::chrono::duration<double> relaxation_time {parameters.get_parameter<double>("relaxation_time")};
    auto simulator = EngineType(restart_file_name);
    auto Calendering_output = simulator.get_output("output_0");
    simulator.remove_output(Calendering_output);

    auto calendering_surface = simulator.get_surface<EngineType::PointSurfacePointer>("top_plate");
    //Remove calendering surface from RVE
    calendering_surface->move(Vec3(0, 0,  1), Vec3(0, 0, 0));

    auto mat = dynamic_cast<ElectrodeMaterial *>(simulator.get_material(0));
//    mat->new_binder_contacts = false;

//==================SET TIME STEP AND MASS SCALING=======================================================================
    double time_step = parameters.get_parameter<double>("time_step")*1E-6;
    std::chrono::duration<double> time_step_us {time_step};
    std::cout << "Time step is:" << time_step_us.count()*1E6 << " µs\n";
    simulator.set_time_increment(time_step_us);
    float mass_scaling_factor = parameters.get_parameter<float>("mass_scaling");
    std::cout << "Mass scaling is:" << mass_scaling_factor << "\n";
    simulator.set_mass_scale_factor(mass_scaling_factor);
//======================================================================================================================

//========================STRAIN RATES==================================================================================
    double strain_rate = parameters.get_parameter<double>("strain_rate");
    std::cout << "Strain rate: " << strain_rate << "\n";
    double strain_point_relaxation = parameters.get_parameter<double>("strain")
//    double strain_point_relaxation = 1.65e-2;
//======================================================================================================================

//====================MAKE OUTPUT PRESCRIBED==============================================================================
    auto output_frequency = parameters.get_parameter<double>("output_frequency");
    std::chrono::duration<double> output_interval {1.0 / output_frequency};
    std::cout << "Output interval is: "<<  output_interval.count() <<"s\n";
    std::cout << "Number of outputs are: "<< relaxation_time.count() / output_interval.count() <<"s\n";
    auto relaxation_output = simulator.create_output(output_directory, output_interval);
    relaxation_output->print_particles = true;
    relaxation_output->print_surface_positions = true;
    relaxation_output->print_kinetic_energy = true;
    relaxation_output->print_contacts = true;
    relaxation_output->print_surface_forces = true;
    relaxation_output->print_fabric_force_tensor =true;
    relaxation_output->print_periodic_bc = true;
    relaxation_output->print_mirror_particles= true;
//======================================================================================================================

    auto deformable_surface = simulator.get_surface<EngineType::DeformablePointSurfacePointer>("bottom_plate");
    /*
    std::cout << "**************** Load step 0 - pre relaxation ****************\n";
    std::chrono::duration<double> pre_relaxation_time {10};
    EngineType::RunForTime run_for_pre_relaxation_time(simulator,pre_relaxation_time);
    simulator.run(run_for_pre_relaxation_time);
    */

    std::cout << "**************** Load step 1 ****************\n";
    std::chrono::duration<double> strain_time_0_relax {(strain_point_relaxation)/strain_rate};
    std::cout << "Loading time: " << strain_time_0_relax.count() << "s\n";
    EngineType::RunForTime run_for_time_BC_stretch(simulator,strain_time_0_relax);
    simulator.set_periodic_boundary_condition_strain_rate('x',loading_direction*strain_rate);
    deformable_surface->set_in_plane_strain_rates(loading_direction*strain_rate,0);
    simulator.run(run_for_time_BC_stretch);

    std::cout << "**************** Hold 1 ****************\n";
    std::cout << "Relaxation time: " << relaxation_time.count() << "s\n";
    simulator.set_periodic_boundary_condition_strain_rate('x',0);
    deformable_surface->set_in_plane_strain_rates(0,0);
    run_for_time_BC_stretch.reset(relaxation_time);
    simulator.run(run_for_time_BC_stretch);
    std::cout << "**************** Simulation finished ****************\n";
}