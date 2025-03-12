//
// Adapted from fracturing_electrode_cycling.cpp by Axel on 2025-02-26
//

#include "../../simulations.h"

#include <vector>

#include "../../../engine/engine.h"
#include "../../../contact_models/Positive_electrode/fracturing_swelling_elastic_plastic_binder_elastic_plastic_particle.h"
#include "../../../materials/electrode_material.h"


void DEM::fracturing_electrode_uniaxial_loading(const std::string &settings_file_name)
{
    using namespace DEM;
    using ForceModel = fracturing_swelling_elastic_plastic_binder_elastic_plastic_particle;
    using ParticleType = FracturableSwellingSphericalParticle<ForceModel>;
    using EngineType = Engine<ForceModel, ParticleType>;
    using namespace std::chrono_literals;

    SimulationParameters parameters(settings_file_name);
    auto restart_file_name = parameters.get_parameter<std::string>("restart_file_name");
    auto output_directory = parameters.get_parameter<std::string>("output_dir");
    auto simulator = EngineType(restart_file_name);
    auto Calendering_output = simulator.get_output("output_0");
    simulator.remove_output(Calendering_output);

    auto mat = dynamic_cast<ElectrodeMaterial *>(simulator.get_material(0));
    mat->new_binder_contacts = false;

    //=REMOVE CALENDERING SURFACE FROM RVE==============================================================================
    auto calendering_surface = simulator.get_surface<EngineType::PointSurfacePointer>("top_plate");
    auto bbox = simulator.get_bounding_box(); //get the XYZ max/min that contain all particles
    double h_1 = bbox[5]; //height of uppermost particle (Z-max)
    double surface_removal_distance = h_1 + 1 - calendering_surface->get_points()[0].z();

    calendering_surface->move(Vec3(0, 0,  surface_removal_distance), Vec3(0, 0, 0));
    //==================================================================================================================

    auto strain_rate = parameters.get_parameter<double>("strain_rate");
    std::cout << "strain_rate =" << strain_rate << "\n";
    auto strain_level = parameters.get_parameter<double>("strain_level");
    std::cout << "strain_level = " << strain_level << "\n";

    if (strain_level < 0.)
    {
        std::cout << "Compression strain applied\n";
        strain_rate = -strain_rate;
    }
    std::chrono::duration<double>loading_time{strain_level/strain_rate};
    std::cout << "loading_time = " << loading_time.count() << "\n";

//==================SET TIME STEP AND MASS SCALING======================================================================
    double time_step = parameters.get_parameter<double>("time_step") * 1e-6;
    std::chrono::duration<double> time_step_us{time_step};
    std::cout << "Time step is:" << time_step_us.count() * 1E6 << " µs\n";
    simulator.set_time_increment(time_step_us);
    float mass_scaling_factor = parameters.get_parameter<float>("mass_scaling");
    std::cout << "Mass scaling is:" << mass_scaling_factor << "\n";
    simulator.set_mass_scale_factor(mass_scaling_factor);
//======================================================================================================================

//====================MAKE OUTPUT PRESCRIBED============================================================================
    double output_number = parameters.get_parameter<double>("output_number");
    std::chrono::duration<double> output_interval{loading_time.count() / output_number};
    std::cout << "Output interval is: " << output_interval.count() << "s\n";
    std::cout << "Number of outputs are: " << loading_time.count() / output_interval.count() << "s\n";
    auto uniaxial_output = simulator.create_output(output_directory, output_interval);
    uniaxial_output->print_particles = true;
    uniaxial_output->print_fractured_particles = true;
    uniaxial_output->print_surface_positions = true;
    uniaxial_output->print_kinetic_energy = true;
    uniaxial_output->print_contacts = true;
    uniaxial_output->print_surface_forces = true;
    uniaxial_output->print_fabric_force_tensor = true;
    uniaxial_output->print_periodic_bc = true;
    uniaxial_output->print_mirror_particles = true;
    //======================================================================================================================

    //==REST LAYER BEFORE SIMULATION====================================================================================
    std::chrono::duration<double> resting_time{loading_time.count() / 10.};
    std::cout << "Running for resting time: " << resting_time.count() << std::endl;
    EngineType::RunForTime run_for_resting_time(simulator, resting_time);
    simulator.run(run_for_resting_time);
    //==================================================================================================================

    //=APPLY BOUNDARY CONDITION LOADING AND RUN SIMULATION==============================================================
    auto deformable_surface = simulator.get_surface<EngineType::DeformablePointSurfacePointer>("bottom_plate");
    deformable_surface->set_in_plane_strain_rates(strain_rate,0.);
    simulator.set_periodic_boundary_condition_strain_rate('x', strain_rate);

    EngineType::RunForTime run_for_loading(simulator, loading_time);
    simulator.run(run_for_loading);
    deformable_surface->set_in_plane_strain_rates(0.,0.);
    simulator.set_periodic_boundary_condition_strain_rate('x', 0.);

    std::cout << "Writing restart file\n";
//    std::stringstream restart_file_name;
//    restart_file_name << "/uniaxial_loaded.res";
//    simulator.write_restart_file(output_directory + restart_file_name.str());
    simulator.write_restart_file(output_directory + "/uniaxial_loaded.res");
    std::cout << "**************** Simulation finalized ****************\n";
    //==================================================================================================================

}