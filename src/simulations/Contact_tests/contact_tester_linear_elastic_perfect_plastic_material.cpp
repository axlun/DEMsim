//
// Created by erolsson on 2019-04-13.
//

#include <chrono>

#include "../simulations.h"

#include "../../engine/engine.h"
#include "../../particles/spherical_particle.h"
#include "../../engine/contact.h"
#include "../../contact_models/stone_material_contact.h"
#include "../../materials/stone_material.h"
#include "../../utilities/file_reading_functions.h"
#include "../../contact_models/viscoelastic.h"
#include "../../contact_models/porous_electrode_contact.h"
#include "../../contact_models/viscoelastic_binder_with_deformable_particles.h"
#include "../../contact_models/elastic_perfect_plastic.h"
#include "../../materials/electrode_material.h"
#include "../../materials/porous_electrode_material.h"



void DEM::contact_tester_elastic_perfect_plastic_material(const std::string& settings_file_name)
{
    using namespace DEM;
    using ForceModel = Elastic_perfectly_plastic_binder;
    using ParticleType = SphericalParticle<ForceModel>;
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
    mat.bt =parameters.get_parameter<double>("bt");
    mat.binder_radius_fraction=parameters.get_parameter<double>("binder_radius_fraction");
    mat.binder_yield_stress_=parameters.get_parameter<double>("binder_yield_stress_");
    std::cout << "Binder thickness, bt:" << mat.bt << std::endl;
    //mat.unloading_exponent = parameters.get_parameter<double>("unloading_exponent");
    mat.fraction_binder_contacts =parameters.get_parameter<double>("fraction_binder_contacts");
    auto p1 = SphericalParticle<ForceModel>(radius, Vec3{-radius-(mat.bt)/2-tick ,0 , 0},
                                            Vec3{}, &mat, 1);
    auto p2 = SphericalParticle<ForceModel>(radius,Vec3{radius+(mat.bt)/2+tick ,0, 0},
                                            Vec3{}, &mat, 1);
    auto c = Contact<ForceModel, ParticleType>(&p2, &p1, timestep);

    //p1.move(Vec3{h1, 0, 0});
    //p2.move(Vec3{-h1, 0, 0});

    fs::path path_to_output_file {filename};
    fs::create_directories(path_to_output_file.parent_path());
    std::ofstream output_file;
    output_file.open(filename);
    auto simulation_time = 0.;
    for(unsigned i = 0; i != increments; ++i) {
        p1.move(Vec3{tick,0 , 0});
        p2.move(Vec3{-tick,0 , 0});
        c.update();
        simulation_time += simulation_time_step;
        output_file << c.get_overlap() << ", " << c.get_normal_force().x() << ", "
                    << p2.get_position().x() - p1.get_position().x() << ", "
                    //  << c.get_tangential_force().y() << ", "  //not any tangential force yet
                    << p1.get_position().y() - p2.get_position().y() << "," << simulation_time << "," << i << std::endl;
    }
    for(unsigned i = 0; i != 2*increments; ++i) {
        p1.move(Vec3{-tick,0 , 0});
        p2.move(Vec3{tick,0 , 0});
        c.update();
        simulation_time += simulation_time_step;
        output_file << c.get_overlap() << ", " << c.get_normal_force().x() << ", "
                    << p2.get_position().x() - p1.get_position().x() << ", "
                    //  << c.get_tangential_force().y() << ", "  //not any tangential force yet
                    << p1.get_position().y() - p2.get_position().y() << "," << simulation_time << "," << i << std::endl;
    }
    for(unsigned i = 0; i != 10*increments; ++i) {
        p1.move(Vec3{tick,0 , 0});
        p2.move(Vec3{-tick,0 , 0});
        c.update();
        simulation_time += simulation_time_step;
        output_file << c.get_overlap() << ", " << c.get_normal_force().x() << ", "
                    << p2.get_position().x() - p1.get_position().x() << ", "
                    //  << c.get_tangential_force().y() << ", "  //not any tangential force yet
                    << p1.get_position().y() - p2.get_position().y() << "," << simulation_time << "," << i << std::endl;
    }
    for(unsigned i = 0; i != 10*increments; ++i) {
        p1.move(Vec3{-tick,0 , 0});
        p2.move(Vec3{tick,0 , 0});
        c.update();
        simulation_time += simulation_time_step;
        output_file << c.get_overlap() << ", " << c.get_normal_force().x() << ", "
                   << p2.get_position().x() - p1.get_position().x() << ", "
                   //  << c.get_tangential_force().y() << ", "  //not any tangential force yet
                   << p1.get_position().y() - p2.get_position().y() << "," << simulation_time << "," << i << std::endl;
    }
}
