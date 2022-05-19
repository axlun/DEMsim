//
// Created by Axel on 2022-05-18.
//

#include "../simulations.h"

#include <vector>

#include "../../engine/engine.h"
#include "../../contact_models/viscoelastic.h"
#include "../../contact_models/Binder_behavour_investigation/viscoelastic_binder_El_Pl_particles.h"
#include "../../materials/electrode_material.h"


void DEM::restart_electrode_calendering(const std::string &settings_file_name)
{
    using namespace DEM;
    using ForceModel = viscoelastic_binder_El_Pl_particles;
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

//    for (auto& p: particles_)
//    {
//    p->radii //get the larges particle radii
//    }
//    simulator.get_particles()

//    for (const auto& p: particles_) {
//        std::cout << p->get_output_string() << "\n";
//    }
    double surface_velocity =parameters.get_parameter<double>("calendaring_surface_velocity"); //How was this chosen?

//=====================================================WRTIE RESTART PROGRAM HERE=======================================
    //Calendaring process
//    std::cout << "****************Initialize calendaring process ****************\n";
//    bbox = simulator.get_bounding_box(); //get the XYZ max/min that contain all particles
//    double h_3 = bbox[5]; //height of uppermost particle (Z-max)
//    std::cout<<"Heigh of uppermoast particle: "<< h_3<< std::endl;
//    if (h_3<1.2*mat->active_particle_height){
//        h_3=1.2*mat->active_particle_height;
//        std::cout<<"Heigh of uppermoast particle lower then 1.2 h_al: "<< h_3<< std::endl;
//    }
//    top_surface->move(-Vec3(0, 0,  mat->active_particle_height*2 - h_3-1.01*max_binder_thickness), Vec3(0, 0, 0)); //Move top surface to uppermost particle+binder thickness
//    top_surface->set_velocity(Vec3(0,0,0.-surface_velocity));
//    std::chrono::duration<double> compaction_time_2 {((h_3+1.01*max_binder_thickness-mat->active_particle_height) / surface_velocity)};
//
//    std::cout<<"Compaction time: "<< compaction_time_2.count()<< std::endl;
//    run_for_time.reset(compaction_time_2);
//    simulator.run(run_for_time);
//    std::cout<<"Writing restart file ";
//    simulator.write_restart_file(output_directory + "/compact_restart_file.res");
//    std::cout << "****************Initialize unloading ****************\n";
//    top_surface->set_velocity(Vec3(0, 0, surface_velocity));
////    mat->adhesive = true; //Enable adhesive when particles are compacted and before unloaded
//    EngineType::SurfaceNormalForceLess zero_force(top_surface, 0.);
//    // simulator.set_rotation(false);
//    simulator.run(zero_force);
//
//    top_surface->set_velocity(Vec3(0, 0, 0));
//    run_for_time.reset(.5s);
//    simulator.run(run_for_time);
//    bbox = simulator.get_bounding_box(); //get the XYZ max/min that contain all particles
//    double Active_layer_height = bbox[5];
//    std::cout<<"Height of active layer "<< Active_layer_height<< std::endl;
//    double NMC_porosity = 1-particle_volume/(box_side*box_side*Active_layer_height) ;
//    std::cout<<"NMC Porosity: "<< NMC_porosity<< std::endl;
//    std::cout<<"Writing restart file ";
//    simulator.write_restart_file(output_directory + "/calendered_electrode_restart_file.res");
//    std::ofstream results_file;
//    results_file.open(output_directory + "/calendered_electrode_results_file.dou");
//    results_file << "Volume of particles=" << particle_volume << "\n";
//    results_file<< "RVE side length= " << box_side << "\n";
//    results_file << "NMC Porosity=" << NMC_porosity << "\n";
//    results_file << "Height of active layer=" << Active_layer_height << "\n";
//    results_file.close();
    }