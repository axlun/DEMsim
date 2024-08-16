//
// Created by Axel on 2024-08-14.
//

#include "../../engine/engine.h"
#include "../../particles/swelling_spherical_particle.h"
#include "../../contact_models/Positive_electrode/swelling_elastic_plastic_binder_elastic_plastic_particle.h"
#include "../../materials/electrode_material.h"

void DEM::material_scaling_particle_tester(const std::string& settings_file_name)
{
    using namespace DEM;
    using ForceModel = swelling_elastic_plastic_binder_elastic_plastic_particle;
    using ParticleType = SwellingSphericalParticle<ForceModel>;
    using EngineType = Engine<ForceModel, ParticleType>;
    SimulationParameters parameters{settings_file_name};

    double mass_scaling = parameters.get_parameter<double>("mass_scaling");
    double time_step = parameters.get_parameter<double>("time_step")*1E-6;
    std::chrono::duration<float> time_step_us {time_step};
    EngineType simulator(time_step_us);

    auto radius_1 = parameters.get_parameter<double>("R1");
    auto radius_2 = parameters.get_parameter<double>("R2");
    auto output_directory = parameters.get_parameter<std::string>("output_dir");
    auto output_number = parameters.get_parameter<double>("output_number");
    auto scaling_rate = parameters.get_parameter<double>("scaling_rate");
    auto normalised_overlap = parameters.get_parameter<double>("normalised_overlap");
    std::cout << "scaling_rate = " << scaling_rate << std::endl;
    std::chrono::duration<double> scaling_time{parameters.get_parameter<double>("scaling_time")};
    std::cout << "scaling time = " << scaling_time.count() << std::endl;

    auto mat = simulator.create_material<ElectrodeMaterial>(4800);
    mat->E = parameters.get_parameter<double>("E");
    mat->nu = parameters.get_parameter<double>("nu");
    mat->Ep=parameters.get_parameter<double>("Ep");
    mat->particle_yield_stress_ = parameters.get_parameter<double>("particle_yield_stress_");
    mat->nup=parameters.get_parameter<double>("nup");

    auto master_curve_file = parameters.get_parameter<std::string>("master_curve_file");
    auto master_curve_vec = read_vector_from_file<double>(master_curve_file);
    mat->F_1_ = master_curve_vec[0];
    mat->alpha_1_ = master_curve_vec[1];
    mat->F_2_ = master_curve_vec[2];
    mat->alpha_2_ = master_curve_vec[3];
    mat->a_1_ = master_curve_vec[4];
    mat->beta_1_ = master_curve_vec[5];
    mat->a_2_ = master_curve_vec[6];
    mat->beta_2_ = master_curve_vec[7];

    mat->binder_radius_fraction=parameters.get_parameter<double>("binder_radius_fraction");
    mat->binder_thickness_fraction=parameters.get_parameter<double>("binder_thickness_fraction");
    mat->binder_yield_stress_=parameters.get_parameter<double>("binder_yield_stress_");
    mat->fraction_binder_contacts =parameters.get_parameter<double>("fraction_binder_contacts");
    mat->alpha_i = parameters.get_vector<double>("alpha_i");
    mat->tau_i =parameters.get_vector<double>("tau_i");
    mat->adhesive = true;

    // ==TIME CONTROL===================================================================================================
//    std::chrono::duration<double> scaling_time(0.2);
    std::chrono::duration<double> move_time = scaling_time * 1E-1;

    double R0_ = 1 / (1 / radius_1 + 1 / radius_2);
    std::cout << "R0: " << R0_ << std::endl;
    auto particle_normal_displacement = R0_ * normalised_overlap;
    auto particle_velocity = particle_normal_displacement/2.0/move_time.count();

    // ==CREATE PARTICLE================================================================================================
    auto p1 = simulator.create_particle(
            radius_1, Vec3{-radius_1*1.001,0,0}, Vec3{particle_velocity,0,0}, mat);
//    p1->set_material_scale_rate(scaling_rate);
    auto p2 = simulator.create_particle(
            radius_2, Vec3{(radius_2)*1.001,0,0}, Vec3{-particle_velocity,0,0}, mat);
//    p2->set_material_scale_rate(scaling_rate);
    simulator.set_mass_scale_factor(mass_scaling);
    double bb_stretch = 1E-1;
    simulator.setup(bb_stretch);

    std::chrono::duration<double>output_duration(scaling_time.count()/output_number);
    std::cout << "output duration is: " << output_duration.count() << "s" << std::endl;
    auto output = simulator.create_output(output_directory,output_duration);
    output->print_particles = true;
    output->print_contacts = true;
    output->print_fabric_force_tensor = true;
    output->print_kinetic_energy = true;

    EngineType::RunForTime RunForTime(simulator,move_time);
    std::cout << "Running for:" << move_time.count() << "s" << std::endl;
    simulator.run(RunForTime);

    for (auto& p: simulator.get_particles())
    {
        p->set_material_scale_rate(scaling_rate);
        p->set_velocity(Vec3(0, 0, 0));
    }

    RunForTime.reset(scaling_time);
    std::cout << "Running for:" << scaling_time.count() << "s" << std::endl;
    simulator.run(RunForTime);

    for (auto& p: simulator.get_particles())
    {
        p->set_material_scale_rate(0);
    }

    p1->set_velocity(Vec3(-particle_velocity, 0, 0));
    p2->set_velocity(Vec3(particle_velocity, 0, 0));
    RunForTime.reset(move_time/2);
    std::cout << "Running for:" << move_time.count() << "s" << std::endl;
    simulator.run(RunForTime);


    for (auto& p: simulator.get_particles())
    {
        p->set_material_scale_rate(-scaling_rate);
        p->set_velocity(Vec3(0, 0, 0));
    }

    RunForTime.reset(scaling_time);
    std::cout << "Running for:" << scaling_time.count() << "s" << std::endl;
    simulator.run(RunForTime);

    for (auto& p: simulator.get_particles())
    {
        p->set_material_scale_rate(0);
    }

    p1->set_velocity(Vec3(particle_velocity, 0, 0));
    p2->set_velocity(Vec3(-particle_velocity, 0, 0));
    RunForTime.reset(move_time);
    std::cout << "Running for:" << move_time.count() << "s" << std::endl;
    simulator.run(RunForTime);

    return;
}

