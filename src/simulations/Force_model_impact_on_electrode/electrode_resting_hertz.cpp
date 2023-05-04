//
// Created by Axel on 2022-05-18.
//

#include "../simulations.h"

#include <vector>

#include "../../engine/engine.h"
#include "../../contact_models/viscoelastic.h"
#include "../../contact_models/Positive_electrode/elastic_plastic_binder_hertz_plastic_particle.h"
#include "../../materials/electrode_material.h"


void DEM::electrode_resting_hertz(const std::string &settings_file_name)
{
    using namespace DEM;
    using ForceModel = elastic_plastic_binder_hertz_plastic_particle;
    using ParticleType = SphericalParticle<ForceModel>;
    using EngineType = Engine<ForceModel,ParticleType>;
    using namespace std::chrono_literals;

    SimulationParameters parameters(settings_file_name);
    auto restart_file_name =parameters.get_parameter<std::string>("restart_file_name");
    auto simulator = EngineType(restart_file_name);
    auto output_directory = parameters.get_parameter<std::string>("output_dir");
    auto mat = simulator.get_material(0);

    auto Old_output = simulator.get_output("output_0");
    simulator.remove_output(Old_output);

    auto calendering_surface = simulator.get_surface<EngineType::PointSurfacePointer>("top_plate");
    //Remove calendering surface from RVE
//    calendering_surface->move(Vec3(0, 0,  1), Vec3(0, 0, 0));


//==================SET TIMESTEP AND MASS SCALING=======================================================================
    double time_step = parameters.get_parameter<double>("time_step")*1e-6;
    std::chrono::duration<double> time_step_us {time_step};
    std::cout << "Time step is:" << time_step_us.count()*1E6 << " µs\n";
    simulator.set_time_increment(time_step_us);
    float mass_scaling_factor = parameters.get_parameter<float>("mass_scaling");
    std::cout << "Mass scaling is:" << mass_scaling_factor << "\n";
    simulator.set_mass_scale_factor(mass_scaling_factor);
//======================================================================================================================

//====================RESTING TIME======================================================================================
    double resting_time = parameters.get_parameter<double>("resting_time");
    std::chrono::duration<double> resting_time_s {resting_time};
    std::cout << "Resting time is: " << resting_time_s.count() << "s \n";
//======================================================================================================================

//====================MAKE OUTPUT DYNAMIC==============================================================================
    std::chrono::duration<double> output_interval {resting_time_s.count()/100};
    std::cout << "Output interval is: "<<  output_interval.count() <<"s\n";
    auto resting_output = simulator.create_output(output_directory, output_interval);
    resting_output->print_particles = true;
    resting_output->print_surface_positions = true;
    resting_output->print_kinetic_energy = true;
    resting_output->print_contacts = true;
    resting_output->print_surface_forces = true;
    resting_output->print_fabric_force_tensor =true;
    resting_output->print_periodic_bc = true;
    resting_output->print_mirror_particles= true;
//======================================================================================================================

//==============RUNNING RESTING TIME====================================================================================
    EngineType::RunForTime run_for_resting_time(simulator,resting_time_s);
    simulator.run(run_for_resting_time);
//======================================================================================================================


//========================STRAIN RATES==================================================================================
//    double strain_rate = parameters.get_parameter<double>("strain_rate");
//    std::cout << "Strain rate: " << strain_rate << "\n";
////    double loading_time_input = parameters.get_parameter<double>("loading_time");
////    std::chrono::duration<double> fix_compaction_time (loading_time_input);
////Strain points in %, converting with e-2
//    double strain_point_0  = 0.00;
//    double strain_point_1  = 0.25e-2;
//    double strain_point_2  = 0.5e-2;
//    double strain_point_3  = 1e-2;
//    double strain_point_4  = 1.5e-2;
//    double strain_point_5  = 2e-2;
//    double unload_strain   = 0.1e-2;
//    std::chrono::duration<double> unload_time {unload_strain/strain_rate};
////======================================================================================================================
//
////====================MAKE OUTPUT DYNAMIC==============================================================================
//    std::chrono::duration<double> output_interval {unload_time.count()/10};
//    std::cout << "Output interval is: "<<  output_interval.count() <<"s\n";
//    auto mechanical_loading_output = simulator.create_output(output_directory, output_interval);
//    mechanical_loading_output->print_particles = true;
//    mechanical_loading_output->print_surface_positions = true;
//    mechanical_loading_output->print_kinetic_energy = true;
//    mechanical_loading_output->print_contacts = true;
//    mechanical_loading_output->print_surface_forces = true;
//    mechanical_loading_output->print_fabric_force_tensor =true;
//    mechanical_loading_output->print_periodic_bc = true;
//    mechanical_loading_output->print_mirror_particles= true;
////======================================================================================================================
//
//    auto deformable_surface = simulator.get_surface<EngineType::DeformablePointSurfacePointer>("bottom_plate");
//    std::cout << "**************** Load step 1 ****************\n";
//
//    std::chrono::duration<double> strain_time_0_1 {(strain_point_1-strain_point_0)/strain_rate};
//    std::cout << "Loading time: " << strain_time_0_1.count() << "s\n";
//    EngineType::RunForTime run_for_time_BC_stretch(simulator,strain_time_0_1);
//    simulator.set_periodic_boundary_condition_strain_rate('x',loading_direction*strain_rate);
//    deformable_surface->set_in_plane_strain_rates(loading_direction*strain_rate,0);
//    simulator.run(run_for_time_BC_stretch);
//
//    std::cout << "**************** Unload and reload 1 ****************\n";
//    std::cout << "Unload and reload time: " << unload_time.count() << "s\n";
//    simulator.set_periodic_boundary_condition_strain_rate('x',-loading_direction*strain_rate);
//    deformable_surface->set_in_plane_strain_rates(-loading_direction*strain_rate,0);
//    run_for_time_BC_stretch.reset(unload_time);
//    simulator.run(run_for_time_BC_stretch);
//    simulator.set_periodic_boundary_condition_strain_rate('x',loading_direction*strain_rate);
//    deformable_surface->set_in_plane_strain_rates(loading_direction*strain_rate,0);
//    run_for_time_BC_stretch.reset(unload_time);
//    simulator.run(run_for_time_BC_stretch);
//
//
//    std::cout << "**************** Load step 2 ****************\n";
//
//    std::chrono::duration<double> strain_time_1_2 {(strain_point_2-strain_point_1)/strain_rate};
//    std::cout << "Loading time: " << strain_time_1_2.count() << "s\n";
//    run_for_time_BC_stretch.reset(strain_time_1_2);
//    simulator.set_periodic_boundary_condition_strain_rate('x',loading_direction*strain_rate);
//    deformable_surface->set_in_plane_strain_rates(loading_direction*strain_rate,0);
//    simulator.run(run_for_time_BC_stretch);
//
//    std::cout << "**************** Unload and reload 2 ****************\n";
//    std::cout << "Unload and reload time: " << unload_time.count() << "s\n";
//    simulator.set_periodic_boundary_condition_strain_rate('x',-loading_direction*strain_rate);
//    deformable_surface->set_in_plane_strain_rates(-loading_direction*strain_rate,0);
//    run_for_time_BC_stretch.reset(unload_time);
//    simulator.run(run_for_time_BC_stretch);
//    simulator.set_periodic_boundary_condition_strain_rate('x',loading_direction*strain_rate);
//    deformable_surface->set_in_plane_strain_rates(loading_direction*strain_rate,0);
//    run_for_time_BC_stretch.reset(unload_time);
//    simulator.run(run_for_time_BC_stretch);
//
//
//
//    std::cout << "**************** Load step 3 ****************\n";
//
//    std::chrono::duration<double> strain_time_2_3 {(strain_point_3-strain_point_2)/strain_rate};
//    std::cout << "Loading time: " << strain_time_2_3.count() << "s\n";
//    run_for_time_BC_stretch.reset(strain_time_2_3);
//    simulator.set_periodic_boundary_condition_strain_rate('x',loading_direction*strain_rate);
//    deformable_surface->set_in_plane_strain_rates(loading_direction*strain_rate,0);
//    simulator.run(run_for_time_BC_stretch);
//
//    std::cout << "**************** Unload and reload 3 ****************\n";
//    std::cout << "Unload and reload time: " << unload_time.count() << "s\n";
//    simulator.set_periodic_boundary_condition_strain_rate('x',-loading_direction*strain_rate);
//    deformable_surface->set_in_plane_strain_rates(-loading_direction*strain_rate,0);
//    run_for_time_BC_stretch.reset(unload_time);
//    simulator.run(run_for_time_BC_stretch);
//    simulator.set_periodic_boundary_condition_strain_rate('x',loading_direction*strain_rate);
//    deformable_surface->set_in_plane_strain_rates(loading_direction*strain_rate,0);
//    run_for_time_BC_stretch.reset(unload_time);
//    simulator.run(run_for_time_BC_stretch);
//
//
//    std::cout << "**************** Load step 4 ****************\n";
//
//    std::chrono::duration<double> strain_time_3_4 {(strain_point_4-strain_point_3)/strain_rate};
//    std::cout << "Loading time: " << strain_time_3_4.count() << "s\n";
//    run_for_time_BC_stretch.reset(strain_time_3_4);
//    simulator.set_periodic_boundary_condition_strain_rate('x',loading_direction*strain_rate);
//    deformable_surface->set_in_plane_strain_rates(loading_direction*strain_rate,0);
//    simulator.run(run_for_time_BC_stretch);
//
//    std::cout << "**************** Unload and reload 4 ****************\n";
//    std::cout << "Unload and reload time: " << unload_time.count() << "s\n";
//    simulator.set_periodic_boundary_condition_strain_rate('x',-loading_direction*strain_rate);
//    deformable_surface->set_in_plane_strain_rates(-loading_direction*strain_rate,0);
//    run_for_time_BC_stretch.reset(unload_time);
//    simulator.run(run_for_time_BC_stretch);
//    simulator.set_periodic_boundary_condition_strain_rate('x',loading_direction*strain_rate);
//    deformable_surface->set_in_plane_strain_rates(loading_direction*strain_rate,0);
//    run_for_time_BC_stretch.reset(unload_time);
//    simulator.run(run_for_time_BC_stretch);
//
//
//    std::cout << "**************** Load step 5 ****************\n";
//
//    std::chrono::duration<double> strain_time_4_5 {(strain_point_5-strain_point_4)/strain_rate};
//    std::cout << "Loading time: " << strain_time_4_5.count() << "s\n";
//    run_for_time_BC_stretch.reset(strain_time_4_5);
//    simulator.set_periodic_boundary_condition_strain_rate('x',loading_direction*strain_rate);
//    deformable_surface->set_in_plane_strain_rates(loading_direction*strain_rate,0);
//    simulator.run(run_for_time_BC_stretch);
//
//    std::cout << "**************** Unload and reload 5 ****************\n";
//    std::cout << "Unload and reload time: " << unload_time.count() << "s\n";
//    simulator.set_periodic_boundary_condition_strain_rate('x',-loading_direction*strain_rate);
//    deformable_surface->set_in_plane_strain_rates(-loading_direction*strain_rate,0);
//    run_for_time_BC_stretch.reset(unload_time);
//    simulator.run(run_for_time_BC_stretch);
//    simulator.set_periodic_boundary_condition_strain_rate('x',loading_direction*strain_rate);
//    deformable_surface->set_in_plane_strain_rates(loading_direction*strain_rate,0);
//    run_for_time_BC_stretch.reset(unload_time);
//    simulator.run(run_for_time_BC_stretch);
//
//
//    std::cout << "**************** Final Unload ****************\n";
//    std::chrono::duration<double> strain_time_5_0 {(strain_point_5-strain_point_0)/strain_rate};
//    std::cout << "Loading time: " << strain_time_5_0.count() << "s\n";
//    run_for_time_BC_stretch.reset(strain_time_5_0);
//    simulator.set_periodic_boundary_condition_strain_rate('x',-loading_direction*strain_rate);
//    deformable_surface->set_in_plane_strain_rates(-loading_direction*strain_rate,0);
//    simulator.run(run_for_time_BC_stretch);
//    std::cout << "**************** Simulation finished ****************\n";
//=====================================================STRETCH THE PERIODIC BCs=======================================
//    std::cout << "**************** Load step 1 ****************\n";
//    EngineType::RunForTime run_for_time_BC_stretch(simulator,.25s);
//    simulator.set_periodic_boundary_condition_strain_rate('x',loading_direction*0.01);
//    deformable_surface->set_in_plane_strain_rates(loading_direction*0.01,0);
//    simulator.run(run_for_time_BC_stretch);
//
//    std::cout << "**************** Unload and reload 1 ****************\n";
//    simulator.set_periodic_boundary_condition_strain_rate('x',-loading_direction*0.01);
//    deformable_surface->set_in_plane_strain_rates(-loading_direction*0.01,0);
//    run_for_time_BC_stretch.reset(.1s);
//    simulator.run(run_for_time_BC_stretch);
//    simulator.set_periodic_boundary_condition_strain_rate('x',loading_direction*0.01);
//    deformable_surface->set_in_plane_strain_rates(loading_direction*0.01,0);
//    run_for_time_BC_stretch.reset(.1s);
//    simulator.run(run_for_time_BC_stretch);
//
//    std::cout << "**************** Load step 2 ****************\n";
//    run_for_time_BC_stretch.reset(.25s);
//    simulator.set_periodic_boundary_condition_strain_rate('x',loading_direction*0.01);
//    deformable_surface->set_in_plane_strain_rates(loading_direction*0.01,0);
//    simulator.run(run_for_time_BC_stretch);
//
//    std::cout << "**************** Unload and reload 2 ****************\n";
//    simulator.set_periodic_boundary_condition_strain_rate('x',-loading_direction*0.01);
//    deformable_surface->set_in_plane_strain_rates(-loading_direction*0.01,0);
//    run_for_time_BC_stretch.reset(.1s);
//    simulator.run(run_for_time_BC_stretch);
//    simulator.set_periodic_boundary_condition_strain_rate('x',loading_direction*0.01);
//    deformable_surface->set_in_plane_strain_rates(loading_direction*0.01,0);
//    run_for_time_BC_stretch.reset(.1s);
//    simulator.run(run_for_time_BC_stretch);
//
//    std::cout << "**************** Load step 3 ****************\n";
//    run_for_time_BC_stretch.reset(.5s);
//    simulator.set_periodic_boundary_condition_strain_rate('x',loading_direction*0.01);
//    deformable_surface->set_in_plane_strain_rates(loading_direction*0.01,0);
//    simulator.run(run_for_time_BC_stretch);
//
//    std::cout << "**************** Unload and reload 3 ****************\n";
//    simulator.set_periodic_boundary_condition_strain_rate('x',-loading_direction*0.01);
//    deformable_surface->set_in_plane_strain_rates(-loading_direction*0.01,0);
//    run_for_time_BC_stretch.reset(.1s);
//    simulator.run(run_for_time_BC_stretch);
//    simulator.set_periodic_boundary_condition_strain_rate('x',loading_direction*0.01);
//    deformable_surface->set_in_plane_strain_rates(loading_direction*0.01,0);
//    run_for_time_BC_stretch.reset(.1s);
//    simulator.run(run_for_time_BC_stretch);
//
//    std::cout << "**************** Load step 4 ****************\n";
//    run_for_time_BC_stretch.reset(.25s);
//    simulator.set_periodic_boundary_condition_strain_rate('x',loading_direction*0.01);
//    deformable_surface->set_in_plane_strain_rates(loading_direction*0.01,0);
//    simulator.run(run_for_time_BC_stretch);
//
//    std::cout << "**************** Unload and reload 4 ****************\n";
//    simulator.set_periodic_boundary_condition_strain_rate('x',-loading_direction*0.01);
//    deformable_surface->set_in_plane_strain_rates(-loading_direction*0.01,0);
//    run_for_time_BC_stretch.reset(.1s);
//    simulator.run(run_for_time_BC_stretch);
//    simulator.set_periodic_boundary_condition_strain_rate('x',loading_direction*0.01);
//    deformable_surface->set_in_plane_strain_rates(loading_direction*0.01,0);
//    run_for_time_BC_stretch.reset(.1s);
//    simulator.run(run_for_time_BC_stretch);
//
//    std::cout << "**************** Load step 5 ****************\n";
//    run_for_time_BC_stretch.reset(.75s);
//    simulator.set_periodic_boundary_condition_strain_rate('x',loading_direction*0.01);
//    deformable_surface->set_in_plane_strain_rates(loading_direction*0.01,0);
//    simulator.run(run_for_time_BC_stretch);
//
//    std::cout << "**************** Unload and reload 5 ****************\n";
//    simulator.set_periodic_boundary_condition_strain_rate('x',-loading_direction*0.01);
//    deformable_surface->set_in_plane_strain_rates(-loading_direction*0.01,0);
//    run_for_time_BC_stretch.reset(.1s);
//    simulator.run(run_for_time_BC_stretch);
//    simulator.set_periodic_boundary_condition_strain_rate('x',loading_direction*0.01);
//    deformable_surface->set_in_plane_strain_rates(loading_direction*0.01,0);
//    run_for_time_BC_stretch.reset(.1s);
//    simulator.run(run_for_time_BC_stretch);
//
//    std::cout << "**************** Final Unload ****************\n";
//    run_for_time_BC_stretch.reset(2s);
//    simulator.set_periodic_boundary_condition_strain_rate('x',-loading_direction*0.01);
//    deformable_surface->set_in_plane_strain_rates(-loading_direction*0.01,0);
//    simulator.run(run_for_time_BC_stretch);



//    std::cout << "**************** resting for 2s ****************\n";
//    EngineType::RunForTime run_for_time(simulator, 2E0s);//Let particles relax for 2s
//    run_for_time_BC_stretch.reset(2s);
//    simulator.run(run_for_time_BC_stretch);
//
//    std::cout << "****************Unloading periodic BCs ****************\n";
//    simulator.set_periodic_boundary_condition_strain_rate('x',0.005);
//    deformable_surface->set_in_plane_strain_rates(0.005,0);
//    run_for_time_BC_stretch.reset(2s);
//    simulator.run(run_for_time_BC_stretch);
//    std::cout << "**************** resting for 2s ****************\n";
//    run_for_time_BC_stretch.reset(2s);
//    simulator.run(run_for_time_BC_stretch);
    }