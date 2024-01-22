//
// Created by Axel on 2023-11-15.
//


#include <chrono>

#include "../../simulations.h"

#include "../../../engine/engine.h"
#include "../../../particles/spherical_particle.h"
#include "../../../surfaces/deformable_point_surface.h"
#include "../../../engine/contact.h"
#include "../../../utilities/file_reading_functions.h"
#include "../../../contact_models/Positive_electrode/elastic_plastic_binder_hertz_plastic_particle.h"
#include "../../../materials/electrode_material.h"
#include "../../../materials/porous_electrode_material.h"

void DEM::particle_force_overlap_hertz_particles(const std::string &settings_file_name)
{
    using namespace DEM;
    using ForceModel = elastic_plastic_binder_hertz_plastic_particle;
    using ParticleType = SphericalParticle<ForceModel>;
    using EngineType = Engine<ForceModel, ParticleType>;
    SimulationParameters parameters{settings_file_name};

    auto mass_scaling = parameters.get_parameter<double>("mass_scaling");
    float time_step = parameters.get_parameter<double>("time_step")*1e-6; //Time input in µs
    std::chrono::duration<double> time_step_us {time_step};
    EngineType simulator(time_step_us); //orig  1E0

    auto mat = simulator.create_material<ElectrodeMaterial>(4800);
    auto radius_1 = parameters.get_parameter<double>("R1");
    auto radius_2 = parameters.get_parameter<double>("R2");
    std::cout << "R1, R2: "<< radius_1 << ", "<< radius_2 << "m" << std::endl;
    auto normalised_overlap = parameters.get_parameter<double>("normalised_overlap");
    auto output_directory = parameters.get_parameter<std::string>("output_dir");
    double output_number = parameters.get_parameter<double>("output_number");
    mat->E = parameters.get_parameter<double>("E");
    mat->nu = parameters.get_parameter<double>("nu");
    mat->Ep=parameters.get_parameter<double>("Ep");
    mat->nup=parameters.get_parameter<double>("nup");
    mat->yield_displacement_coeff=parameters.get_parameter<double>("yield_displacement_coeff");

    mat->binder_radius_fraction=parameters.get_parameter<double>("binder_radius_fraction");
    mat->binder_thickness_fraction=parameters.get_parameter<double>("binder_thickness_fraction");
    mat->binder_yield_stress_=parameters.get_parameter<double>("binder_yield_stress_");
    mat->fraction_binder_contacts =parameters.get_parameter<double>("fraction_binder_contacts");
    mat->alpha_i = parameters.get_vector<double>("alpha_i");
    mat->tau_i =parameters.get_vector<double>("tau_i");
    mat->adhesive = true;

// =PARTICLE CONTACT ===================================================================================================

    auto p1 = simulator.create_particle(radius_1,Vec3{-radius_1*1.00010,0 , 0},Vec3{0,0,0}, mat);
    auto p2 = simulator.create_particle(radius_2,Vec3{radius_2*1.00010,0 , 0},Vec3{0,0,0}, mat);
    simulator.set_mass_scale_factor(mass_scaling);

    simulator.setup();

//=====Time controlled=================================================================================================================
    std::chrono::duration<double> move_time( 1.0);
    double R0_ = 1 / (1 / radius_1 + 1 / radius_2);
    std::cout << "R0: "<< R0_<< "m\n";

    auto particle_normal_displacement =  R0_ * normalised_overlap;
    auto particle_velocity = particle_normal_displacement/2.0/move_time.count();

    std::chrono::duration<double> output_duration(move_time.count()/output_number);

    auto contact_output = simulator.create_output(output_directory,output_duration);
    contact_output->print_particles = true;
    contact_output->print_contacts = true;
    contact_output->print_fabric_force_tensor = true;
    contact_output->print_kinetic_energy = true;

    std::cout << "Particle normal disp: "<< particle_normal_displacement<< "m\n";
    std::cout << "Move_time is: "<< move_time.count()<< "s\n";
    std::cout << "Output_duration is: "<< output_duration.count()<< "s\n";

    EngineType ::RunForTime RunForTime(simulator,move_time);
//    simulator.run(RunForTime);

//Normal movement
    double num_array[] = {1.0, -1.0};//{.2, -.2, .9, -.8, .9};

    for (double n : num_array)
    {
        p1->set_velocity(Vec3{particle_velocity*n, 0,0});
        p2->set_velocity(Vec3{-particle_velocity*n, 0, 0});
        RunForTime.reset(move_time);
        simulator.run(RunForTime);
    }
}
// =END OF SCRIPTi======================================================================================================
