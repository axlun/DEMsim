//
// Created by Axel on 2021-11-19.
//

#include <chrono>

#include "../simulations.h"

#include "../../engine/engine.h"
#include "../../particles/spherical_particle.h"
#include "../../surfaces/deformable_point_surface.h"
#include "../../engine/contact.h"
#include "../../utilities/file_reading_functions.h"
#include "../../contact_models/viscoelastic.h"
#include "../../contact_models/porous_electrode_contact.h"
#include "../../contact_models/viscoelastic_binder_with_deformable_particles.h"
#include "../../contact_models/elastic_perfect_plastic.h"
//#include "../../contact_models/Binder_behavour_investigation/elastic_plastic_binder_rigid_perfect_plastic_particle.h"
#include "../../contact_models/Binder_behavour_investigation/elastic_plastic_binder_hertz_plastic_particle.h"
#include "../../materials/electrode_material.h"
#include "../../materials/porous_electrode_material.h"

void DEM::particle_contact_tester_viscoelastic_binder_El_Pl_particles(const std::string &settings_file_name)
{
    using namespace DEM;
    using ForceModel = elastic_plastic_binder_hertz_plastic_particle;
    //    using ForceModel = elastic_plastic_binder_rigid_perfect_plastic_particle;
    using ParticleType = SphericalParticle<ForceModel>;
    using SurfaceType = DeformablePointSurface<ForceModel,ParticleType>;
    using namespace std::chrono_literals;
    namespace fs = std::filesystem;

    using Material = ElectrodeMaterial;
    Material mat {0, 4800.};

    SimulationParameters parameters{settings_file_name};

    auto radius = parameters.get_parameter<double>("R");
    auto increments = parameters.get_parameter<unsigned>("N");
    auto simulation_time_step = parameters.get_parameter<double>("t");
    std::chrono::duration<double> timestep = std::chrono::duration<double>{parameters.get_parameter<double>("t")};
    //auto h1 = parameters.get_parameter<double>("h1");
    std::cout << "time step:" << timestep.count() <<"s"<< std::endl;
    auto tick = parameters.get_parameter<double>("tick");
    std::cout << "tick:" << tick << std::endl;
    auto filename= parameters.get_parameter<std::string>("output_file");
    mat.E = parameters.get_parameter<double>("E");
    mat.nu = parameters.get_parameter<double>("nu");
    mat.Ep=parameters.get_parameter<double>("Ep");
    mat.nup=parameters.get_parameter<double>("nup");
    //mat.bt =parameters.get_parameter<double>("bt");
    mat.binder_radius_fraction=parameters.get_parameter<double>("binder_radius_fraction");
    mat.binder_thickness_fraction=parameters.get_parameter<double>("binder_thickness_fraction");
    mat.binder_stiffness_coefficient=parameters.get_parameter<double>("binder_stiffness_coefficient");
    mat.yield_displacement_coeff=parameters.get_parameter<double>("yield_displacement_coeff");
    mat.binder_yield_stress_=parameters.get_parameter<double>("binder_yield_stress_");
    mat.particle_yield_stress_=parameters.get_parameter<double>("particle_yield_stress_");

    //std::cout << "Binder thickness, bt:" << mat.bt << std::endl;
    //mat.unloading_exponent = parameters.get_parameter<double>("unloading_exponent");
    mat.fraction_binder_contacts =parameters.get_parameter<double>("fraction_binder_contacts");
    mat.alpha_i = parameters.get_vector<double>("alpha_i");
    mat.tau_i =parameters.get_vector<double>("tau_i");
   // For binder contact
//    auto p1 = SphericalParticle<ForceModel>(radius, Vec3{-radius-(mat.binder_thickness_fraction*radius)/2-4*tick ,0 , 0},
//                                            Vec3{}, &mat, 1);
//    auto p2 = SphericalParticle<ForceModel>(radius,Vec3{radius+(mat.binder_thickness_fraction*radius)/2+4*tick ,0, 0},
//                                            Vec3{}, &mat, 1);
//  For particle contact
    auto p1 = SphericalParticle<ForceModel>(radius, Vec3{-radius-1*tick ,0 , 0},
                                            Vec3{}, &mat, 1);
    auto p2 = SphericalParticle<ForceModel>(radius,Vec3{radius+1*tick ,0, 0},
                                            Vec3{}, &mat, 1);
    auto c = Contact<ForceModel, ParticleType>(&p2, &p1, timestep);

    p1.move(Vec3{tick, 0, 0});
    p2.move(Vec3{-tick, 0, 0});

    fs::path path_to_output_file {filename};
    fs::create_directories(path_to_output_file.parent_path());
    std::ofstream output_file;
    output_file.open(filename);
    auto simulation_time = 0.;
    //Loading of particle
    for(unsigned i = 0; i != increments; ++i) {
        p1.move(Vec3{tick,0 , 0});
        p2.move(Vec3{-tick,0 , 0});
        c.update();
        simulation_time += simulation_time_step;
//        std::cout << "Force:" << c.get_normal_force().x()<<"N"<< std::endl;
//        std::cout << "P1 pos:" << p1.get_position().x()<<""<< std::endl;
//        std::cout << "P2 pos:" << p2.get_position().x()<<""<< std::endl;

        output_file << c.get_overlap() << ", " << c.get_normal_force().x() << ", "
                    << p2.get_position().x() - p1.get_position().x() << ", "
                    //  << c.get_tangential_force().y() << ", "  //not any tangential force yet
                    << p1.get_position().y() - p2.get_position().y() << "," << simulation_time << ", " << i << ", " << c.get_normal() << std::endl;
    }
    //Unloading of particle
    for(unsigned i = 0; i != .2*increments; ++i) {
        p1.move(Vec3{-tick,0 , 0});
        p2.move(Vec3{tick,0 , 0});
        c.update();
        simulation_time += simulation_time_step;
        output_file << c.get_overlap() << ", " << c.get_normal_force().x() << ", "
                    << p2.get_position().x() - p1.get_position().x() << ", "
                    //  << c.get_tangential_force().y() << ", "  //not any tangential force yet
                    << p1.get_position().y() - p2.get_position().y() << "," << simulation_time << ", " << i << ", " << c.get_normal() << std::endl;
    }
    //Reload of particle contact
    for(unsigned i = 0; i != .4*increments; ++i) {
        p1.move(Vec3{tick,0 , 0});
        p2.move(Vec3{-tick,0 , 0});
        c.update();
        simulation_time += simulation_time_step;
        output_file << c.get_overlap() << ", " << c.get_normal_force().x() << ", "
                    << p2.get_position().x() - p1.get_position().x() << ", "
                    //  << c.get_tangential_force().y() << ", "  //not any tangential force yet
                    << p1.get_position().y() - p2.get_position().y() << "," << simulation_time << ", " << i << ", " << c.get_normal() << std::endl;
    }
    //Unloading of particle contact
    for(unsigned i = 0; i != .5*increments; ++i) {
        p1.move(Vec3{-tick,0 , 0});
        p2.move(Vec3{tick,0 , 0});
        c.update();
        simulation_time += simulation_time_step;
        output_file << c.get_overlap() << ", " << c.get_normal_force().x() << ", "
                    << p2.get_position().x() - p1.get_position().x() << ", "
                    //  << c.get_tangential_force().y() << ", "  //not any tangential force yet
                    << p1.get_position().y() - p2.get_position().y() << "," << simulation_time << ", " << i << ", " << c.get_normal() << std::endl;
    }
//    //Reloading of particle contact
//    for(unsigned i = 0; i != 1*increments; ++i) {
//        p1.move(Vec3{tick,0 , 0});
//        p2.move(Vec3{-tick,0 , 0});
//        c.update();
//        simulation_time += simulation_time_step;
//        output_file << c.get_overlap() << ", " << c.get_normal_force().x() << ", "
//                    << p2.get_position().x() - p1.get_position().x() << ", "
//                    //  << c.get_tangential_force().y() << ", "  //not any tangential force yet
//                    << p1.get_position().y() - p2.get_position().y() << "," << simulation_time << ", " << i << ", " << c.get_normal() << std::endl;
//    }
//    //Unloading of particle contact
//    for(unsigned i = 0; i != 1*increments; ++i) {
//        p1.move(Vec3{-tick,0 , 0});
//        p2.move(Vec3{tick,0 , 0});
//        c.update();
//        simulation_time += simulation_time_step;
//        output_file << c.get_overlap() << ", " << c.get_normal_force().x() << ", "
//                    << p2.get_position().x() - p1.get_position().x() << ", "
//                    //  << c.get_tangential_force().y() << ", "  //not any tangential force yet
//                    << p1.get_position().y() - p2.get_position().y() << "," << simulation_time << ", " << i << ", " << c.get_normal() << std::endl;
//    }

    //    for(unsigned i = 0; i != 1*increments; ++i) {
//        p1.move(Vec3{tick,0 , 0});
//        p2.move(Vec3{-tick,0 , 0});
//        c.update();
//        simulation_time += simulation_time_step;
//    output_file << c.get_overlap() << ", " << c.get_normal_force().x() << ", "
//                << p2.get_position().x() - p1.get_position().x() << ", "
//                //  << c.get_tangential_force().y() << ", "  //not any tangential force yet
//                << p1.get_position().y() - p2.get_position().y() << "," << simulation_time << ", " << i << ", " << c.get_normal() << std::endl;
//    }
//    for(unsigned i = 0; i != 1*increments; ++i) {
//        p1.move(Vec3{-tick,0 , 0});
//        p2.move(Vec3{tick,0 , 0});
//        c.update();
//        simulation_time += simulation_time_step;
//    output_file << c.get_overlap() << ", " << c.get_normal_force().x() << ", "
//                << p2.get_position().x() - p1.get_position().x() << ", "
//                //  << c.get_tangential_force().y() << ", "  //not any tangential force yet
//                << p1.get_position().y() - p2.get_position().y() << "," << simulation_time << ", " << i << ", " << c.get_normal() << std::endl;
//    }


//==================================================//
//                Surface contact                   //
//==================================================//
//    auto surf_point_1 = Vec3(0, 1, 1);
//    auto surf_point_2 = Vec3( 0, -1, 1);
//    auto surf_point_3 = Vec3( 0,  -1, -1);
//    auto surf_point_4 = Vec3(0,  1, -1);
//    std::vector<Vec3> surface_points{surf_point_1, surf_point_4, surf_point_2, surf_point_3};
//    auto deformable_surface = DeformablePointSurface<ForceModel, ParticleType>( 1, surface_points, true, "Def_Surface",
//                                                                true, 0);
//
////  auto deformable_surface = simulator.create_deformable_point_surface(surface_points, true);
//
//    auto p3 = SphericalParticle<ForceModel>(radius, Vec3{+radius+(mat.binder_thickness_fraction*radius)+4*tick ,0 , 0},
//                                            Vec3{}, &mat, 1);
//
//    auto surf_contact = Contact<ForceModel, ParticleType>(&p3, &deformable_surface, timestep);
//
//    for(unsigned i = 0; i != 6*increments; ++i) {
//        p3.move(Vec3{-tick,0 , 0});
//        surf_contact.update(); //comment
//        simulation_time += simulation_time_step;
//        output_file << surf_contact.get_overlap() << ", " << surf_contact.get_normal_force().x() << ", "
//                    << p3.get_position().x() << ", "
//                    //  << c.get_tangential_force().y() << ", "  //not any tangential force yet
//                    << p3.get_position().y() << "," << simulation_time << "," << i << "," << surf_contact.get_normal() << std::endl;
//    }
//    for(unsigned i = 0; i != 10*increments; ++i) {
//        //p3.move(Vec3{-tick,0 , 0});
//        surf_contact.update(); //comment
//        simulation_time += simulation_time_step;
//        output_file << surf_contact.get_overlap() << ", " << surf_contact.get_normal_force().x() << ", "
//                    << p3.get_position().x() << ", "
//                    //  << c.get_tangential_force().y() << ", "  //not any tangential force yet
//                    << p3.get_position().y() << "," << simulation_time << "," << i << "," << surf_contact.get_normal() << std::endl;
//    }
//    for(unsigned i = 0; i != 8*increments; ++i) {
//        p3.move(Vec3{+tick,0 , 0});
//        surf_contact.update(); //comment
//        simulation_time += simulation_time_step;
//        output_file << surf_contact.get_overlap() << ", " << surf_contact.get_normal_force().x() << ", "
//                    << p3.get_position().x() << ", "
//                    //  << c.get_tangential_force().y() << ", "  //not any tangential force yet
//                    << p3.get_position().y() << "," << simulation_time << "," << i << "," << surf_contact.get_normal() << std::endl;
//    }
//    for(unsigned i = 0; i != 10*increments; ++i) {
//        //p3.move(Vec3{-tick,0 , 0});
//        surf_contact.update(); //comment
//        simulation_time += simulation_time_step;
//        output_file << surf_contact.get_overlap() << ", " << surf_contact.get_normal_force().x() << ", "
//                    << p3.get_position().x() << ", "
//                    //  << c.get_tangential_force().y() << ", "  //not any tangential force yet
//                    << p3.get_position().y() << "," << simulation_time << "," << i << "," << surf_contact.get_normal() << std::endl;
//    }
}