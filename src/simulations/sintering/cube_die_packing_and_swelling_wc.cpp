//
// Adapted from swelling_cube_die.cpp on 2024-05-28
//
// Created by Axel on 2024-03-21
//

#include <string>

#include "../../engine/engine.h"
#include "../../materials/electrode_material.h"
#include "../../surfaces/point_surface.h"
#include "../../utilities/file_reading_functions.h"
#include "../../utilities/filling_functions.h"
#include "../simulations.h"
#include "../../contact_models/Positive_electrode/swelling_elastic_plastic_binder_elastic_plastic_particle.h"

void DEM::cube_die_packing_and_swelling_wc(const std::string& settings_file_name)
{
    using namespace DEM;
    using ForceModel = swelling_elastic_plastic_binder_elastic_plastic_particle;
    using ParticleType = SwellingSphericalParticle<ForceModel>;
    using SurfaceType = DeformablePointSurface<ForceModel, ParticleType>;
    using EngineType = Engine<ForceModel, ParticleType>;
    SimulationParameters parameters{settings_file_name};

    double time_step = parameters.get_parameter<double>("time_step") * 1E-6;
    std::chrono::duration<double> time_step_us{time_step};
    std::cout << "Time step is: " << time_step_us.count() << " s" << std::endl;
    EngineType simulator(time_step_us);

    double mass_scaling = parameters.get_parameter<double>("mass_scaling");
    simulator.set_mass_scale_factor(mass_scaling);
    std::cout << "Mass scaling is: " << mass_scaling << std::endl;

    auto output_directory = parameters.get_parameter<std::string>("output_dir");
    auto output_number = parameters.get_parameter<double>("output_number");

    auto delta_t = parameters.get_parameter<double>("delta_t");

    auto N = parameters.get_parameter<std::size_t>("N");
    auto wc_radius = parameters.get_parameter<double>("wc_diameter")/2;

    auto wc_mat = simulator.create_material<ElectrodeMaterial>(15630);
    auto wc_alpha = parameters.get_parameter<double>("wc_alpha");
    auto wc_master_curve_file = parameters.get_parameter<std::string>("wc_master_curve_file");
    auto wc_master_curve_vec = read_vector_from_file<double>(wc_master_curve_file);
    wc_mat->F_1_ = wc_master_curve_vec[0];
    wc_mat->alpha_1_ = wc_master_curve_vec[1];
    wc_mat->F_2_ = wc_master_curve_vec[2];
    wc_mat->alpha_2_ = wc_master_curve_vec[3];
    wc_mat->a_1_ = wc_master_curve_vec[4];
    wc_mat->beta_1_ = wc_master_curve_vec[5];
    wc_mat->a_2_ = wc_master_curve_vec[6];
    wc_mat->beta_2_ = wc_master_curve_vec[7];
    wc_mat->Ep = parameters.get_parameter<double>("wc_Ep");
    wc_mat->nup = parameters.get_parameter<double>("wc_nup");
    wc_mat->rhop = parameters.get_parameter<double>("wc_rhop");
    wc_mat->particle_yield_stress_ = parameters.get_parameter<double>("wc_particle_yield_stress_");
    auto alpha_2 = parameters.get_parameter<double>("alpha_2");

    auto particle_density_at_filling = parameters.get_parameter<double>("filling_density");
    auto pre_packing_density= parameters.get_parameter<double>("pre_packing_density");
    auto packing_density= parameters.get_parameter<double>("packing_density");


    std::chrono::duration<double> pre_packing_time {parameters.get_parameter<double>("pre_packing_time")};
    std::chrono::duration<double> packing_time {parameters.get_parameter<double>("packing_time")};
    std::chrono::duration<double> resting_time {parameters.get_parameter<double>("resting_time")};
    std::chrono::duration<double> swelling_time {parameters.get_parameter<double>("swelling_time")};
    std::chrono::duration<double> resting_time_2 {parameters.get_parameter<double>("resting_time_2")};
    auto wc_swell_rate = wc_alpha * delta_t * wc_radius / swelling_time.count();
    auto swell_rate_2 = alpha_2 * delta_t * wc_radius / swelling_time.count();

    // ================================================================================================================
    // Creating the particles
    // ================================================================================================================
    std::vector<double> particle_radii (N, wc_radius);

    double max_radii = *std::max_element(particle_radii.begin(), particle_radii.end());
    double particle_volume = 0.;
    for( const auto &r: particle_radii)
    {
        particle_volume += 4./3. * 3.1415 * r * r * r;
    }

    std::cout << "Volume of particles are: " << particle_volume << std::endl;
    //============================================================================================================
    //Creates the cube so that the volume corresponds to an initial density of filling_density

    auto box_side = pow(particle_volume/particle_density_at_filling, 1./3);
    std::cout << "box_side " <<box_side << "\n";

    // ================================================================================================================
    //Creating the box
    // ================================================================================================================

    //Creates all the points that define the initial box
    auto p1 = Vec3(-box_side/2, -box_side/2, -box_side/2);
    auto p2 = Vec3( box_side/2, -box_side/2, -box_side/2);
    auto p3 = Vec3( box_side/2,  box_side/2, -box_side/2);
    auto p4 = Vec3(-box_side/2,  box_side/2, -box_side/2);
    auto p5 = Vec3(-box_side/2, -box_side/2,  box_side/2);
    auto p6 = Vec3( box_side/2, -box_side/2,  box_side/2);
    auto p7 = Vec3( box_side/2,  box_side/2,  box_side/2);
    auto p8 = Vec3(-box_side/2,  box_side/2,  box_side/2);

    //Saves the points corresponding to a surface in a vector
    std::vector<Vec3> front_points{p5, p6, p2, p1};
    std::vector<Vec3> back_points{p7, p8, p4, p3};
    std::vector<Vec3> right_points{p6, p7, p3, p2};
    std::vector<Vec3> left_points{p8, p5, p1, p4};
    std::vector<Vec3> bottom_points{p1, p2, p3, p4};
    std::vector<Vec3> top_points{p8, p7, p6, p5};

    // Creates all surfaces as point surfaces and writes out their normal-vector in the terminal
    std::cout << "Generating surfaces"<<"\n";
    auto top_surface = simulator.create_point_surface(
            top_points, true, "top_surface", false);
    std::cout << "Normal of top surface: " << top_surface->get_normal() << "\n";
    auto bottom_surface = simulator.create_point_surface(
            bottom_points, true, "bottom_surface", false);
    std::cout << "Normal of bottom surface: " << bottom_surface->get_normal() << "\n";
    auto front_surface = simulator.create_point_surface(
            front_points, true, "front_surface", false);
    std::cout << "Normal of front surface: " << front_surface->get_normal() << "\n";
    auto back_surface = simulator.create_point_surface(
            back_points, true, "back_surface", false);
    std::cout << "Normal of back surface: " << back_surface->get_normal() << "\n";
    auto right_surface = simulator.create_point_surface(
            right_points, true, "right_surface", false);
    std::cout << "Normal of right surface: " << right_surface->get_normal() << "\n";
    auto left_surface = simulator.create_point_surface(
            left_points, true, "left_surface", false);
    std::cout << "Normal of left surface: " << left_surface->get_normal() << "\n";

    std::cout << "Surfaces generated"<<"\n";

    //************************************************
    auto particle_positions = random_fill_box(-box_side/2, box_side/2, -box_side/2, box_side/2,
                                              -box_side/2, box_side/2, particle_radii);
    std::cout << "Particle positions generated"<<"\n";
    //************************************************

    //Creates particles at particle positions
    for (std::size_t i = 0; i != particle_positions.size(); ++i)
    {
            simulator.create_particle(wc_radius, particle_positions[i], Vec3(0,0,0), wc_mat);
    }

    std::chrono::duration<double> output_time { swelling_time/output_number};
    auto output1 = simulator.create_output(output_directory , output_time, "output");
    output1->print_particles = true;
    output1->print_kinetic_energy = true;
    output1->print_surface_positions = true;
    output1->print_surface_forces = true;
    output1->print_fabric_force_tensor = true;
    output1->print_contacts = true;

    //========================
    //Running the Simulation
    //========================
    simulator.setup();
    EngineType::RunForTime run_for_time(simulator, pre_packing_time);
    //====================================
    // Step 0: compress to pre packing density
    //====================================
    double h_target = pow(particle_volume/pre_packing_density, 1./3)/2;
    double surface_velocity = (h_target - top_surface->get_points()[0].z())/(pre_packing_time.count());
    std::cout << "Top surface position is : " << top_surface->get_points()[0].z() << "\n";
    std::cout << "Plate target is: " << h_target << "\n";
    std::cout << "Surface velocity is: " << surface_velocity << "\n";
    //===============Right Direction of surface disp===============
    top_surface->set_velocity(Vec3(0, 0, surface_velocity));
    bottom_surface->set_velocity(Vec3(0, 0, -surface_velocity));
    front_surface->set_velocity(Vec3(0, -surface_velocity, 0));
    back_surface->set_velocity(Vec3(0, surface_velocity, 0));
    left_surface->set_velocity(Vec3(-surface_velocity, 0, 0));
    right_surface->set_velocity(Vec3(surface_velocity, 0, 0));
    //=============================================================

    std::cout << "\nStep 0: Particle pre compaction, running for t = " << pre_packing_time.count() << std::endl;
    simulator.run(run_for_time);

    //====================================
    // Step 1: Compress to pre packing density
    //====================================
    h_target = pow(particle_volume/packing_density, 1./3)/2;
    surface_velocity = (h_target - top_surface->get_points()[0].z())/(packing_time.count());
    std::cout << "Top surface position is : " << top_surface->get_points()[0].z() << "\n";
    std::cout << "Plate target is: " << h_target << "\n";
    std::cout << "Surface velocity is: " << surface_velocity << "\n";
    //===============Right Direction of surface disp===============
    top_surface->set_velocity(Vec3(0, 0, surface_velocity));
    bottom_surface->set_velocity(Vec3(0, 0, -surface_velocity));
    front_surface->set_velocity(Vec3(0, -surface_velocity, 0));
    back_surface->set_velocity(Vec3(0, surface_velocity, 0));
    left_surface->set_velocity(Vec3(-surface_velocity, 0, 0));
    right_surface->set_velocity(Vec3(surface_velocity, 0, 0));
    //=============================================================
    std::cout << "\nStep 1: Particle compaction, running for t = " << packing_time.count() << std::endl;
    run_for_time.reset(packing_time);
    simulator.run(run_for_time);


    //====================================
    // Step 2: Stop particles and rest
    //====================================
    std::cout << "\nSimulation step 2: particle resting, running for t = " << resting_time.count() << std::endl;
    //===============Right Direction of surface disp===============
    top_surface->set_velocity(Vec3(0, 0, 0));
    bottom_surface->set_velocity(Vec3(0, 0, 0));
    front_surface->set_velocity(Vec3(0, 0, 0));
    back_surface->set_velocity(Vec3(0, 0, 0));
    left_surface->set_velocity(Vec3(0, 0, 0));
    right_surface->set_velocity(Vec3(0, 0, 0));
    //=============================================================
    for (auto& p:simulator.get_particles())
    {
        p->set_velocity(Vec3(0,0,0));
    }
    run_for_time.reset(resting_time);
    simulator.run(run_for_time);
    //========================
    // Step 3: swell particles
    //========================
    std::cout << "\nSimulation step 3: particle swelling, running for t = " << swelling_time.count() << std::endl;
    for (auto& p:simulator.get_particles())
    {
        p->set_velocity(Vec3(0,0,0));
        p->set_swell_rate(wc_swell_rate);
        if (p->get_id() % 2 == 0 && alpha_2 > 1e-30)
        {
            std::cout << p->get_id() << std::endl;
            p->set_swell_rate(swell_rate_2);
        }
    }
    run_for_time.reset(swelling_time);
    simulator.run(run_for_time);
    //========================
    // Step 4: rest particles
    //========================
    for (auto& p:simulator.get_particles())
    {
        p->set_swell_rate(0.);
    }
    run_for_time.reset(resting_time_2);
    simulator.run(run_for_time);
    std::cout << "Simulation ended!" << std::endl;
}