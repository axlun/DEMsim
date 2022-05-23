//
// Created by Axel on 2022-05-18.
//

#include "../simulations.h"

#include <vector>

#include "../../engine/engine.h"
#include "../../contact_models/viscoelastic.h"
#include "../../contact_models/Binder_behavour_investigation/viscoelastic_binder_El_Pl_particles.h"
#include "../../materials/electrode_material.h"


void DEM::electrode_mechanical_loading(const std::string &settings_file_name)
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
//    auto mat = simulator.make_material_from_restart_data()
//    auto mat = simulator.create_material<ElectrodeMaterial>(4800);//**************OLD MATERIAL SHOULD BE READ FROM RESTART******
    auto Calendering_output = simulator.get_output("output_0");
    simulator.remove_output(Calendering_output);
    auto mechanical_loading_output = simulator.create_output(output_directory, 0.05s);
    mechanical_loading_output->print_particles = true;
    mechanical_loading_output->print_surface_positions = true;
    mechanical_loading_output->print_kinetic_energy = true;
    mechanical_loading_output->print_contacts = true;
    mechanical_loading_output->print_surface_forces = true;
    mechanical_loading_output->print_fabric_force_tensor =true;
    mechanical_loading_output->print_periodic_bc = true;
    mechanical_loading_output->print_mirror_particles= true;


    auto deformable_surface = simulator.get_surface<EngineType::DeformablePointSurfacePointer>("bottom_plate");


    auto mat = simulator.get_material(0);
    std::cout << "material density" << mat->density << "\n";

//=====================================================STRETCH THE PERIODIC BCs=======================================
    simulator.set_mass_scale_factor(1E2);
    EngineType::RunForTime run_for_time_BC_stretch(simulator,1s);

    simulator.set_periodic_boundary_condition_strain_rate('x',-0.01);
    deformable_surface->set_in_plane_strain_rates(-0.01,0);
    simulator.run(run_for_time_BC_stretch);



//    mat->binder_thickness_fraction = parameters.get_parameter<double>("binder_thickness_fraction");
//    std::cout << "material id0 " << mat->binder_thickness_fraction << "\n";


    //Get the largest binder thickness
//    auto max_binder_thickness = 0.;
//    for (auto& p: simulator.get_particles())
//    {
//        std::cout << "particle radius: " << p->get_radius() << "\n";
//        std::cout << "binder thickness fraction: " << mat->binder_thickness_fraction << "\n";
//        if (mat->binder_thickness_fraction*p->get_radius() > max_binder_thickness)
//        {
//            max_binder_thickness = mat->binder_thickness_fraction*p->get_radius();
//        }//get the larges particle radii
//    }
//    auto deformable_surface = simulator.get_surface<EngineType::DeformablePointSurfacePointer>("bottom_plate");
//
//
//
//
//    std::cout << "****************Time to start the simulation**************** \n";
//
//
////    double surface_velocity =parameters.get_parameter<double>("calendaring_surface_velocity"); //How was this chosen?
//    std::cout << "max_binder_thickness: " << max_binder_thickness << "\n";
//
//


    }