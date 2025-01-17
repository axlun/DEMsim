//
// Created by Axel on 2025-01-13
// based on material_scaling_periodic_compaction.cpp
//

#include "../../engine/engine.h"
#include "../../materials/electrode_material.h"
#include "../../contact_models/Positive_electrode/fracturing_swelling_elastic_plastic_binder_elastic_plastic_particle.h"

void DEM::fracturing_particle_plate_compression(const std::string& settings_file_name)
{
    using namespace DEM;
    using ForceModel = fracturing_swelling_elastic_plastic_binder_elastic_plastic_particle;
    using ParticleType = FracturableSwellingSphericalParticle<ForceModel>;
    using SurfaceType = DeformablePointSurface<ForceModel, ParticleType>;
    using EngineType = Engine<ForceModel, ParticleType>;
    SimulationParameters parameters{settings_file_name};

    double time_step = parameters.get_parameter<double>("time_step") * 1E-6;
    std::chrono::duration<double> time_step_us{time_step};
    std::cout << "Time step is: " << time_step_us.count() << " µs" << std::endl;
    EngineType simulator(time_step_us);

    double mass_scaling = parameters.get_parameter<double>("mass_scaling");
    simulator.set_mass_scale_factor(mass_scaling);
    std::cout << "Mass scaling is: " << mass_scaling << std::endl;

    auto output_directory = parameters.get_parameter<std::string>("output_dir");
    auto output_number = parameters.get_parameter<double>("output_number");

    auto particle_radius = parameters.get_parameter<double>("R1");

    auto mat = simulator.create_material<ElectrodeMaterial>(4800);
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
    mat->E = parameters.get_parameter<double>("E");
    mat->Ep = parameters.get_parameter<double>("Ep");
    mat->nu = parameters.get_parameter<double>("nu");
    mat->nup = parameters.get_parameter<double>("nup");
    mat->rhop = parameters.get_parameter<double>("rhop");
    mat->tau_i = parameters.get_vector<double>("tau_i");
    mat->alpha_i = parameters.get_vector<double>("alpha_i");
    mat->particle_yield_stress_ = parameters.get_parameter<double>("particle_yield_stress_");
    mat->binder_radius_fraction = parameters.get_parameter<double>("binder_radius_fraction");
    mat->binder_thickness_fraction = parameters.get_parameter<double>("binder_thickness_fraction");
    mat->binder_yield_stress_ = parameters.get_parameter<double>("binder_yield_stress_");
    mat->fraction_binder_contacts = parameters.get_parameter<double>("fraction_binder_contacts");

    mat->fracture_degradation_factor_ = parameters.get_parameter<double>("fracture_degradation_factor_");
    mat->particle_fracture_strength_ = parameters.get_parameter<double>("particle_fracture_strength_");

    auto normalised_compression_overlap = parameters.get_parameter<double>("normalised_compression_overlap");

    std::chrono::duration<double> simulation_time{parameters.get_parameter<double>("simulation_time")};
    // ================================================================================================================
    // Creating the particle
    // ================================================================================================================
    auto particle = simulator.create_particle(particle_radius, Vec3(0, 0, particle_radius), Vec3(0, 0, 0), mat);

    // ================================================================================================================
    //Creating the surfaces
    // ================================================================================================================
    //Creates all the points that define the initial box
    auto p1 = Vec3(-particle_radius, -particle_radius, 0);
    auto p2 = Vec3(particle_radius, -particle_radius, 0);
    auto p3 = Vec3(particle_radius, particle_radius, 0);
    auto p4 = Vec3(-particle_radius, particle_radius, 0);
    auto p5 = Vec3(-particle_radius, -particle_radius, particle_radius * 2);
    auto p6 = Vec3(particle_radius, -particle_radius, particle_radius * 2);
    auto p7 = Vec3(particle_radius, particle_radius, particle_radius * 2);
    auto p8 = Vec3(-particle_radius, particle_radius, particle_radius * 2);

    //Saves the points corresponding to a surface in a vector
    std::vector <Vec3> bottom_points{p1, p2, p3, p4};
    std::vector <Vec3> top_points{p8, p7, p6, p5};

    // Creates all surfaces as point surfaces and writes out their normal-vector in the terminal
    std::cout << "Generating surfaces" << "\n";
    auto top_surface = simulator.create_point_surface(
            top_points, true, "top_surface", false);
    std::cout << "Normal of top surface: " << top_surface->get_normal() << "\n";
    auto bottom_surface = simulator.create_point_surface(
            bottom_points, true, "bottom_surface", false);
    std::cout << "Normal of bottom surface: " << bottom_surface->get_normal() << "\n";


    std::chrono::duration<double> output_time{simulation_time / output_number};
    auto output1 = simulator.create_output(output_directory, output_time, "output");
    output1->print_particles = true;
    output1->print_kinetic_energy = true;
    output1->print_surface_positions = true;
    output1->print_surface_forces = true;
    output1->print_fabric_force_tensor = true;
    output1->print_contacts = true;
    output1->print_periodic_bc = true;
    output1->print_mirror_particles = true;

    //========================
    //Running the Simulation
    //========================
    simulator.setup(2*particle_radius);
    EngineType::RunForTime run_for_time(simulator, simulation_time);

    //Set surface velocity and run for compression time
    double surface_velocity = -2 * particle_radius * normalised_compression_overlap/ simulation_time.count();
    top_surface->set_velocity(Vec3(0, 0, surface_velocity));
    simulator.run(run_for_time);

    top_surface->set_velocity(Vec3(0, 0, -surface_velocity));
    run_for_time.reset(simulation_time);
    simulator.run(run_for_time);
}