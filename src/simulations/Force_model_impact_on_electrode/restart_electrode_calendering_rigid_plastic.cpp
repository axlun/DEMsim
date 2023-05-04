//
// Created by Axel on 2022-03-01.
//

#include "../simulations.h"

#include <vector>

#include "../../engine/engine.h"
#include "../../contact_models/viscoelastic.h"
#include "../../contact_models/Positive_electrode/elastic_plastic_binder_rigid_plastic_particle.h"
#include "../../materials/electrode_material.h"

void DEM::restart_electrode_calendering_rigid_plastic(const std::string &settings_file_name)
{
    using namespace DEM;
    using ForceModel = elastic_plastic_binder_rigid_plastic_particle;
    using ParticleType = SphericalParticle<ForceModel>;
    using EngineType = Engine<ForceModel,ParticleType>;
    using namespace std::chrono_literals;

    SimulationParameters parameters(settings_file_name);

    auto restart_file_name =parameters.get_parameter<std::string>("restart_file_name");
    auto output_directory = parameters.get_parameter<std::string>("output_dir");
    auto calendering_height = parameters.get_parameter<double>("calendering_height");
    auto simulator = EngineType(restart_file_name);
    auto mat = dynamic_cast<const ElectrodeMaterial *>(simulator.get_material(0));

//====================TIME STEP AND MASS SCALING========================================================================
    double time_step = parameters.get_parameter<double>("time_step")*1e-6;
    std::chrono::duration<double> time_step_us {time_step};
    std::cout << "Time step is:" << time_step_us.count()*1E6 << " µs\n";
    simulator.set_time_increment(time_step_us);
    float mass_scaling_factor = parameters.get_parameter<float>("mass_scaling");
    std::cout << "Mass scaling is:" << mass_scaling_factor << "\n";
    simulator.set_mass_scale_factor(mass_scaling_factor); //Orig 1E2
//======================================================================================================================

//=======================REMOVE OLD OUTPUT==============================================================================
    auto Calendering_output = simulator.get_output("output_0");
    simulator.remove_output(Calendering_output);
//======================================================================================================================

    auto top_surface = simulator.get_surface<EngineType::PointSurfacePointer>("top_plate");
    auto deformable_surface = simulator.get_surface<EngineType::DeformablePointSurfacePointer>("bottom_plate");
    std::cout << "Binder thickness fraction: " << mat->binder_thickness_fraction << "\n";
    auto max_binder_thickness = 0.;
    double particle_volume = 0.;
    for (auto& p: simulator.get_particles())
    {
        particle_volume += 4. / 3. * pi * p->get_radius() * p->get_radius() * p->get_radius();
        if (mat->binder_thickness_fraction*p->get_radius() > max_binder_thickness)
        {
            max_binder_thickness = mat->binder_thickness_fraction*p->get_radius();
        }//get the larges particle radii
    }
    auto box_side = abs(2*top_surface->get_points()[0].x());
    auto initial_calender_surface_height = top_surface->get_points()[0].z();
    std::cout << "max_binder_thickness: " << max_binder_thickness << "\n";
    std::cout << "Particle volume: " << particle_volume << "\n";
    std::cout << "Box side: " << box_side << "\n";

    //Calendaring process
    auto bbox = simulator.get_bounding_box(); //get the XYZ max/min that contain all particles
    double h_3 = bbox[5]; //height of uppermost particle (Z-max)
    std::cout<<"Heigh of uppermoast particle: "<< h_3<< std::endl;
    std::cout<<"Calendering height : "<< calendering_height << std::endl;
    if (h_3<1.2*calendering_height){
        std::cout<<"Heigh of uppermoast particle lower then 1.2 h_al: "<< h_3 << std::endl;
        h_3=1.2*calendering_height;
        std::cout<<"Calednering surface changed to: "<< h_3 << std::endl;
    }
    EngineType::RunForTime run_for_time(simulator, .1E-0s);
    std::cout << "****************Initialize calendaring process ****************\n";
    top_surface->move(-Vec3(0, 0,  initial_calender_surface_height - h_3-1.01*max_binder_thickness), Vec3(0, 0, 0)); //Move top surface to uppermost particle+binder thickness

//================================SET COMPACTION TIME AND FIND SURF. VELOCITY===========================================
    double compaction_time_input = parameters.get_parameter<double>("compaction_time");
    std::chrono::duration<double> fix_compaction_time (compaction_time_input);
    auto calendering_distance = (h_3+1.01*max_binder_thickness-calendering_height);
    auto fix_time_surface_vel = calendering_distance/fix_compaction_time.count();
    std::cout << "Calendering surface velocity: " << fix_time_surface_vel << "m/s \n";

//============================Calculate output frequency================================================================
    double output_number =parameters.get_parameter<double>("output_number");
    std::chrono::duration<double> output_interval {fix_compaction_time.count()/output_number};
    std::cout << "Number of outputs: " << output_number << "\n";
    std::cout << "Output interval: " << output_interval.count() << "s\n";

    auto restart_output = simulator.create_output(output_directory, output_interval);
    restart_output->print_particles = true;
    restart_output->print_surface_positions = true;
    restart_output->print_kinetic_energy = true;
    restart_output->print_contacts = true;
    restart_output->print_surface_forces = true;
    restart_output->print_fabric_force_tensor =true;
    restart_output->print_periodic_bc = true;
    restart_output->print_mirror_particles= true;
//======================================================================================================================

    run_for_time.reset(fix_compaction_time);
    top_surface->set_velocity(Vec3(0,0,-fix_time_surface_vel));
    std::cout << "Calendering time: " << fix_compaction_time.count() << "s\n";
    simulator.run(run_for_time);

//======================================================================================================================

//=======================UNLOADING CALENDERING SURFACE TO ZERO FORCE=================================================
    std::cout << "****************Unloading untill zero force on calendering surface ****************\n";
    EngineType::SurfaceNormalForceLess zero_force(top_surface, 0.);
    top_surface->set_velocity(Vec3(0,0,fix_time_surface_vel));
    simulator.run(zero_force);
    top_surface->set_velocity(Vec3(0, 0, 0));
//======================================================================================================================

    auto top_surface_position_after_calendering = top_surface->get_points()[0].z();
    std::cout<<"Top surface position after calendering: "<< top_surface_position_after_calendering<< std::endl;
    std::cout << "****************Let rest for half of calendering time ****************\n";
    std::cout << "Running for: " << (fix_compaction_time*0.5).count() << "s \n";

    run_for_time.reset(fix_compaction_time*0.5);
    simulator.run(run_for_time);
    bbox = simulator.get_bounding_box(); //get the XYZ max/min that contain all particles
    double Active_layer_height = bbox[5];
    std::cout<<"Height of active layer "<< Active_layer_height<< std::endl;
    //Calculate this from the position when the force on calendering surface is 0
    double NMC_porosity = 1-particle_volume/(box_side*box_side*top_surface_position_after_calendering);
//    double NMC_porosity = 1-particle_volume/(box_side*box_side*Active_layer_height) ;
    std::cout<<"NMC Porosity: "<< NMC_porosity<< std::endl;
    std::cout<<"Writing restart file ";
    simulator.write_restart_file(output_directory + "/calendered_electrode_restart_file.res");
    std::ofstream results_file;
    results_file.open(output_directory + "/calendered_electrode_results_file.dou");
    results_file << "Volume of particles=" << particle_volume << "\n";
    results_file<< "RVE side length= " << box_side << "\n";
    results_file << "NMC Porosity=" << NMC_porosity << "\n";
    results_file << "Height of active layer=" << Active_layer_height << "\n";
    results_file << "Top surface position=" << top_surface_position_after_calendering << "\n";
    results_file.close();

    }