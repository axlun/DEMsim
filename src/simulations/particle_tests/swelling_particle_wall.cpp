//
// Created by Axel on 2024-02-15.
//

#include "../../engine/engine.h"
#include "../../particles/swelling_spherical_particle.h"
#include "../../surfaces/deformable_point_surface.h"
#include "../../contact_models/Positive_electrode/swelling_elastic_plastic_binder_elastic_plastic_particle.h"
#include "../../materials/electrode_material.h"

void DEM::swelling_particle_wall_tester(const std::string& settings_file_name)
{
    using namespace DEM;
    using ForceModel = swelling_elastic_plastic_binder_elastic_plastic_particle;
    using ParticleType = SwellingSphericalParticle<ForceModel>;
    using SurfaceType = DeformablePointSurface<ForceModel, ParticleType>;
    using EngineType = Engine<ForceModel, ParticleType>;
    SimulationParameters parameters{settings_file_name};

    double mass_scaling = parameters.get_parameter<double>("mass_scaling");
    double time_step = parameters.get_parameter<double>("time_step")*1E-6;
    std::chrono::duration<float> time_step_us {time_step};
    EngineType simulator(time_step_us);

    auto radius_1 = parameters.get_parameter<double>("R1");
    auto output_directory = parameters.get_parameter<std::string>("output_dir");
    auto output_number = parameters.get_parameter<double>("output_number");
    auto swell_rate = parameters.get_parameter<double>("swell_rate");
    std::cout << "swell_rate = " << swell_rate << std::endl;
    auto numer_of_falls = parameters.get_parameter<double>("number_of_falls");

    auto mat= simulator.create_material<ElectrodeMaterial>(4800);
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
//    mat->binder_stiffness_coefficient=parameters.get_parameter<double>("binder_stiffness_coefficient");
    mat->binder_yield_stress_=parameters.get_parameter<double>("binder_yield_stress_");
    mat->fraction_binder_contacts =parameters.get_parameter<double>("fraction_binder_contacts");
    mat->alpha_i = parameters.get_vector<double>("alpha_i");
    mat->tau_i =parameters.get_vector<double>("tau_i");
    mat->adhesive = true;
    // ==CREATE PARTICLE================================================================================================
    double particle_center_position = radius_1 * 1.101;
    auto p1 = simulator.create_particle(
            radius_1, Vec3{0,0,particle_center_position}, Vec3{0,0,0}, mat);
    p1->set_swell_rate(swell_rate);
//    auto p2 = simulator.create_particle(
//            radius_2, Vec3{(radius_1+radius_2)*1.001,0,0}, Vec3{0,0,0}, mat);
//    p2->set_swell_rate(swell_rate);

    // ==CREATE SURFACE=================================================================================================
    auto Wp1 = Vec3(2*radius_1, 2*radius_1, 0);
    auto Wp2 = Vec3(-2*radius_1, 2*radius_1, 0);
    auto Wp3 = Vec3(-2*radius_1,  -2*radius_1, 0);
    auto Wp4 = Vec3(2*radius_1,  -2*radius_1, 0);
    std::vector<Vec3> wall_points{Wp1, Wp2, Wp3, Wp4};
    auto contact_surface = simulator.create_point_surface(
            wall_points, true, "contact_surface", true);
    std::cout << "Normal of contact surface: " << contact_surface->get_normal() << "\n";

    simulator.set_mass_scale_factor(mass_scaling);
    double bb_stretch = 1E-1;
    simulator.setup(bb_stretch);

    // ==TIME CONTROL===================================================================================================
    double gravity = 9.82;
    simulator.set_gravity(Vec3(0,0, -gravity));
    std::chrono::duration<double> fall_time(pow(2*(particle_center_position-radius_1)/gravity, 0.5));
    std::cout << "Fall time is: " << fall_time.count() << std::endl;
    std::chrono::duration<double> simulation_time( 2 * numer_of_falls * fall_time);

    //std::chrono::duration<double>swell_time(0.2);
    //std::chrono::duration<double>output_duration(swell_time.count()/output_number);
    std::chrono::duration<double>output_duration(simulation_time.count()/output_number);
    auto output = simulator.create_output(output_directory,output_duration);
    output->print_particles = true;
    output->print_contacts = true;
    output->print_fabric_force_tensor = true;
    output->print_kinetic_energy = true;
    output->print_surface_positions = true;
    output->print_surface_forces = true;

    EngineType::RunForTime RunForTime(simulator,simulation_time);
    simulator.run(RunForTime);
    return;
}

