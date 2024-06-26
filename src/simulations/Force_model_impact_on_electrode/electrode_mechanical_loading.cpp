//
// Created by Axel on 2022-05-18.
//

#include "../simulations.h"

#include <vector>

#include "../../engine/engine.h"
#include "../../contact_models/viscoelastic.h"
#include "../../contact_models/Positive_electrode/elastic_plastic_binder_rigid_perfect_plastic_particle_OLD.h"
#include "../../materials/electrode_material.h"


void DEM::electrode_mechanical_loading(const std::string &settings_file_name)
{
    using namespace DEM;
    using ForceModel = elastic_plastic_binder_rigid_perfect_plastic_particle_OLD;
    using ParticleType = SphericalParticle<ForceModel>;
    using EngineType = Engine<ForceModel,ParticleType>;
    using namespace std::chrono_literals;

    SimulationParameters parameters(settings_file_name);
    auto restart_file_name =parameters.get_parameter<std::string>("restart_file_name");
    auto output_directory = parameters.get_parameter<std::string>("output_dir");
    auto loading_direction = parameters.get_parameter<int>("loading_direction"); // compression -> -1, tension -> 1
    auto simulator = EngineType(restart_file_name);
//    auto mat = simulator.make_material_from_restart_data()
//    auto mat = simulator.create_material<ElectrodeMaterial>(4800);//**************OLD MATERIAL SHOULD BE READ FROM RESTART******
    auto Calendering_output = simulator.get_output("output_0");
    simulator.remove_output(Calendering_output);
    auto mechanical_loading_output = simulator.create_output(output_directory, 0.005s); //
    mechanical_loading_output->print_particles = true;
    mechanical_loading_output->print_surface_positions = true;
    mechanical_loading_output->print_kinetic_energy = true;
    mechanical_loading_output->print_contacts = true;
    mechanical_loading_output->print_surface_forces = true;
    mechanical_loading_output->print_fabric_force_tensor =true;
    mechanical_loading_output->print_periodic_bc = true;
    mechanical_loading_output->print_mirror_particles= true;


    auto deformable_surface = simulator.get_surface<EngineType::DeformablePointSurfacePointer>("bottom_plate");
    auto calendering_surface = simulator.get_surface<EngineType::PointSurfacePointer>("top_plate");
    //Remove calendering surface from RVE
    calendering_surface->move(Vec3(0, 0,  1), Vec3(0, 0, 0));

    auto mat = simulator.get_material(0);
//    std::cout << "material density" << mat->density << "\n";

//=====================================================STRETCH THE PERIODIC BCs=======================================
    simulator.set_mass_scale_factor(1E2);

    std::cout << "**************** Load step 1 ****************\n";
    EngineType::RunForTime run_for_time_BC_stretch(simulator,.25s);
    simulator.set_periodic_boundary_condition_strain_rate('x',loading_direction*0.01);
    deformable_surface->set_in_plane_strain_rates(loading_direction*0.01,0);
    simulator.run(run_for_time_BC_stretch);

    std::cout << "**************** Unload and reload 1 ****************\n";
    simulator.set_periodic_boundary_condition_strain_rate('x',-loading_direction*0.01);
    deformable_surface->set_in_plane_strain_rates(-loading_direction*0.01,0);
    run_for_time_BC_stretch.reset(.1s);
    simulator.run(run_for_time_BC_stretch);
    simulator.set_periodic_boundary_condition_strain_rate('x',loading_direction*0.01);
    deformable_surface->set_in_plane_strain_rates(loading_direction*0.01,0);
    run_for_time_BC_stretch.reset(.1s);
    simulator.run(run_for_time_BC_stretch);

    std::cout << "**************** Load step 2 ****************\n";
    run_for_time_BC_stretch.reset(.25s);
    simulator.set_periodic_boundary_condition_strain_rate('x',loading_direction*0.01);
    deformable_surface->set_in_plane_strain_rates(loading_direction*0.01,0);
    simulator.run(run_for_time_BC_stretch);

    std::cout << "**************** Unload and reload 2 ****************\n";
    simulator.set_periodic_boundary_condition_strain_rate('x',-loading_direction*0.01);
    deformable_surface->set_in_plane_strain_rates(-loading_direction*0.01,0);
    run_for_time_BC_stretch.reset(.1s);
    simulator.run(run_for_time_BC_stretch);
    simulator.set_periodic_boundary_condition_strain_rate('x',loading_direction*0.01);
    deformable_surface->set_in_plane_strain_rates(loading_direction*0.01,0);
    run_for_time_BC_stretch.reset(.1s);
    simulator.run(run_for_time_BC_stretch);

    std::cout << "**************** Load step 3 ****************\n";
    run_for_time_BC_stretch.reset(.5s);
    simulator.set_periodic_boundary_condition_strain_rate('x',loading_direction*0.01);
    deformable_surface->set_in_plane_strain_rates(loading_direction*0.01,0);
    simulator.run(run_for_time_BC_stretch);

    std::cout << "**************** Unload and reload 3 ****************\n";
    simulator.set_periodic_boundary_condition_strain_rate('x',-loading_direction*0.01);
    deformable_surface->set_in_plane_strain_rates(-loading_direction*0.01,0);
    run_for_time_BC_stretch.reset(.1s);
    simulator.run(run_for_time_BC_stretch);
    simulator.set_periodic_boundary_condition_strain_rate('x',loading_direction*0.01);
    deformable_surface->set_in_plane_strain_rates(loading_direction*0.01,0);
    run_for_time_BC_stretch.reset(.1s);
    simulator.run(run_for_time_BC_stretch);

    std::cout << "**************** Load step 4 ****************\n";
    run_for_time_BC_stretch.reset(.25s);
    simulator.set_periodic_boundary_condition_strain_rate('x',loading_direction*0.01);
    deformable_surface->set_in_plane_strain_rates(loading_direction*0.01,0);
    simulator.run(run_for_time_BC_stretch);

    std::cout << "**************** Unload and reload 4 ****************\n";
    simulator.set_periodic_boundary_condition_strain_rate('x',-loading_direction*0.01);
    deformable_surface->set_in_plane_strain_rates(-loading_direction*0.01,0);
    run_for_time_BC_stretch.reset(.1s);
    simulator.run(run_for_time_BC_stretch);
    simulator.set_periodic_boundary_condition_strain_rate('x',loading_direction*0.01);
    deformable_surface->set_in_plane_strain_rates(loading_direction*0.01,0);
    run_for_time_BC_stretch.reset(.1s);
    simulator.run(run_for_time_BC_stretch);

    std::cout << "**************** Load step 5 ****************\n";
    run_for_time_BC_stretch.reset(.75s);
    simulator.set_periodic_boundary_condition_strain_rate('x',loading_direction*0.01);
    deformable_surface->set_in_plane_strain_rates(loading_direction*0.01,0);
    simulator.run(run_for_time_BC_stretch);

    std::cout << "**************** Unload and reload 5 ****************\n";
    simulator.set_periodic_boundary_condition_strain_rate('x',-loading_direction*0.01);
    deformable_surface->set_in_plane_strain_rates(-loading_direction*0.01,0);
    run_for_time_BC_stretch.reset(.1s);
    simulator.run(run_for_time_BC_stretch);
    simulator.set_periodic_boundary_condition_strain_rate('x',loading_direction*0.01);
    deformable_surface->set_in_plane_strain_rates(loading_direction*0.01,0);
    run_for_time_BC_stretch.reset(.1s);
    simulator.run(run_for_time_BC_stretch);

    std::cout << "**************** Final Unload ****************\n";
    run_for_time_BC_stretch.reset(2s);
    simulator.set_periodic_boundary_condition_strain_rate('x',-loading_direction*0.01);
    deformable_surface->set_in_plane_strain_rates(-loading_direction*0.01,0);
    simulator.run(run_for_time_BC_stretch);



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