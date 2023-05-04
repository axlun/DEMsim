//
// Created by Axel on 2022-08-30.
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
//#include "../../contact_models/Positive_electrode/elastic_plastic_binder_rigid_perfect_plastic_particle_OLD.h"
#include "../../contact_models/Positive_electrode/elastic_plastic_binder_hertz_plastic_particle.h"
#include "../../materials/electrode_material.h"
#include "../../materials/porous_electrode_material.h"

void DEM::binder_wall_tangential_contact_tester_elastic_plastic_binder_hertz_particles(const std::string &settings_file_name)
{
    using namespace DEM;
    using ForceModel = elastic_plastic_binder_hertz_plastic_particle;
    using ParticleType = SphericalParticle<ForceModel>;
    using EngineType = Engine<ForceModel, ParticleType>;
    using SurfaceType = DeformablePointSurface<ForceModel,ParticleType>;
    using namespace std::chrono_literals;
    namespace fs = std::filesystem;
    SimulationParameters parameters{settings_file_name};


    auto mass_scaling = parameters.get_parameter<double>("mass_scaling");
    float time_step = parameters.get_parameter<double>("time_step")*1e-6;
    std::chrono::duration<double> time_step_us {time_step};
    EngineType simulator(time_step_us); //orig  1E0

    auto mat = simulator.create_material<ElectrodeMaterial>(4800);
//    using Material = ElectrodeMaterial;
//    Material mat {0, 4800.};

    auto radius = parameters.get_parameter<double>("R");
    auto increments = parameters.get_parameter<unsigned>("N");
    auto simulation_time_step = parameters.get_parameter<double>("t");
//    auto particle_velocity = parameters.get_parameter<double>("particle_velocity");
    auto output_directory = parameters.get_parameter<std::string>("output_dir");
    double output_number = parameters.get_parameter<double>("output_number");
    mat->E = parameters.get_parameter<double>("E");
    mat->nu = parameters.get_parameter<double>("nu");
    mat->Ep=parameters.get_parameter<double>("Ep");
    mat->nup=parameters.get_parameter<double>("nup");
    mat->binder_radius_fraction=parameters.get_parameter<double>("binder_radius_fraction");
    mat->binder_thickness_fraction=parameters.get_parameter<double>("binder_thickness_fraction");
    mat->binder_stiffness_coefficient=parameters.get_parameter<double>("binder_stiffness_coefficient");
    mat->yield_displacement_coeff=parameters.get_parameter<double>("yield_displacement_coeff");
    mat->binder_yield_stress_=parameters.get_parameter<double>("binder_yield_stress_");
    mat->particle_yield_stress_=parameters.get_parameter<double>("particle_yield_stress_");

    mat->fraction_binder_contacts =parameters.get_parameter<double>("fraction_binder_contacts");
    mat->alpha_i = parameters.get_vector<double>("alpha_i");
    mat->tau_i =parameters.get_vector<double>("tau_i");
    mat->adhesive = true;

    // ==PARTICLE========================================================================================================
    auto p1 = simulator.create_particle(radius,Vec3{-radius-(mat->binder_thickness_fraction*radius)*1.010,0 , 0},Vec3{0,0,0}, mat);
    auto p2 = simulator.create_particle(radius,Vec3{5*radius+(mat->binder_thickness_fraction*radius)*1.010/2.0,0 , 0},Vec3{0,0,0}, mat);

    // ==WALL===========================================================================================================
    auto Wp1 = Vec3(0, 2*radius, 2*radius);
    auto Wp2 = Vec3(0, 2*radius, -2*radius);
    auto Wp3 = Vec3(0, -2*radius,  -2*radius);
    auto Wp4 = Vec3(0, -2*radius,  2*radius);
    std::vector<Vec3> wall_points{Wp1, Wp2, Wp3, Wp4};


    auto contact_surface = simulator.create_point_surface(wall_points, true,"contact_surface",true);


    simulator.set_mass_scale_factor(mass_scaling);
    simulator.setup(radius*(1+mat->binder_thickness_fraction));
    simulator.set_rotation(false);

/*
    p1->set_velocity(Vec3{particle_velocity,0,0});
    p2->set_velocity(Vec3{-particle_velocity,0,0});
    auto particle_normal_displacement = radius*mat->binder_thickness_fraction*0.01; //Move particles 1% of the binder distance
    std::chrono::duration<double> particle_normal_displacement_time(50*particle_normal_displacement/particle_velocity); //Divide by 2 as both particles are moving
    EngineType ::RunForTime run_for_time(simulator,particle_normal_displacement_time);
    std::cout << "Normal movement of particles: sim time is:" << particle_normal_displacement_time.count() <<"\n";
    simulator.run(run_for_time);
*/

//=====Velocity controlled=================================================================================================================
//    p1->set_velocity(Vec3{0,particle_velocity,0});
//    p2->set_velocity(Vec3{0,-particle_velocity,0});
//    auto particle_tangential_displacement =  radius*mat->binder_thickness_fraction*0.01; //Move particles 1% of the binder distance
//    std::chrono::duration<double>particle_tangential_displacement_time(particle_tangential_displacement/particle_velocity);
//    run_for_time.reset(particle_tangential_displacement_time);
//    std::cout << "Tangential movement of particles"<< "\n";
//    simulator.run(run_for_time);

//=====Time controlled displacement=====================================================================================
    std::chrono::duration<double> move_time( mat->tau_i[0] / 10.0 );
    auto particle_tangential_displacement =  radius*mat->binder_thickness_fraction*0.01; //Move particles 1% of the binder distance
    auto particle_velocity = particle_tangential_displacement/move_time.count();


    // =CREATEING OUTPUT================================================================================================
    std::chrono::duration<double> output_duration(mat->tau_i[0]/output_number);
    auto contact_output = simulator.create_output(output_directory,output_duration);
    contact_output->print_particles = true;
    contact_output->print_contacts = true;
    contact_output->print_fabric_force_tensor = true;
    contact_output->print_kinetic_energy = true;
    contact_output->print_surface_positions = true;
    contact_output->print_surface_forces= true;

    std::cout << "Particle tangential disp: "<< particle_tangential_displacement<< "m\n";
    std::cout << "Move_time is: "<< move_time.count()<< "s\n";
    std::cout << "Output_duration is: "<< output_duration.count()<< "s\n";


    //make particle come into contact
    auto normal_move_distance = (mat->binder_thickness_fraction*radius)*0.010000001;

    auto normal_vel = normal_move_distance/move_time.count();
    p1->set_velocity(Vec3{normal_vel,0,0});
    EngineType ::RunForTime RunForTime(simulator,move_time);
    simulator.run(RunForTime);

    RunForTime.reset(std::chrono::duration<double>( 20*mat->tau_i[0] ));
    simulator.run(RunForTime);


//    normal_move_distance = (mat->binder_thickness_fraction*radius)*0.00000001;
//    normal_vel = normal_move_distance/2.0/move_time.count();
//    std::cout << "normal_move_distance is: "<< normal_move_distance<< "s\n";
//    std::cout << "normal_vel is: "<< normal_vel<< "s\n";
//
//    p1->set_velocity(Vec3{-normal_vel,0,0});
//    p2->set_velocity(Vec3{normal_vel,-0,0});
//    RunForTime.reset(move_time);
//    simulator.run(RunForTime);


//Tangential movement
    double num_array[] = {1};

    for (double n : num_array) {
        p1->set_velocity(Vec3{0,particle_velocity*n,0});
        RunForTime.reset(move_time);
        simulator.run(RunForTime);

//        p1->set_velocity(Vec3{0,0,0});
//        p2->set_velocity(Vec3{0,-0,0});

        RunForTime.reset(std::chrono::duration<double>( 5*mat->tau_i[0] ));
        simulator.run(RunForTime);
    }
//
////    Diagonal movement
    p1->set_velocity(Vec3{0,-0,particle_velocity});
//    p2->set_velocity(Vec3{0,particle_velocity,-particle_velocity});
    RunForTime.reset(move_time);
    simulator.run(RunForTime);
//
//    p1->set_velocity(Vec3{0,0,0});
//    p2->set_velocity(Vec3{0,-0,0});
//
    RunForTime.reset(std::chrono::duration<double>( 40*mat->tau_i[0] ));
    simulator.run(RunForTime);
}