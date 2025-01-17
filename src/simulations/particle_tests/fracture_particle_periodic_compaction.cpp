//
// Created by Axel on 2025-01-15
// based on material_scaling_periodic_compaction.cpp
//

#include "../../engine/engine.h"
#include "../../materials/electrode_material.h"
#include "../../contact_models/Positive_electrode/fracturing_swelling_elastic_plastic_binder_elastic_plastic_particle.h"

void DEM::fracturing_particle_periodic_compaction(const std::string& settings_file_name) {
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

    auto N = parameters.get_parameter<std::size_t>("N");
    auto particle_file = parameters.get_parameter<std::string>("radius_file");

    // Reads a vector with different relative densities where unloading will be made
    auto density_levels = parameters.get_vector<double>( "density_levels" );

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
    // Imports the initial relative density of the particles at filling
    auto particle_density_at_filling = parameters.get_parameter<double>("filling_density");

    std::chrono::duration<double> simulation_time {parameters.get_parameter<double>("simulation_time")};
    // ================================================================================================================
    // Creating the particles
    // ================================================================================================================
    auto particle_radii = read_vector_from_file<double>(particle_file);
    particle_radii.assign(N,*particle_radii.begin());
    if (particle_radii.size() == 1)
    {
        particle_radii.assign(N, *particle_radii.begin());
    }
    else
    {
        particle_radii.assign(particle_radii.begin(), particle_radii.begin() + N);
        std::sort(particle_radii.rbegin(), particle_radii.rend());
    }
    double max_radii = *std::max_element(particle_radii.begin(), particle_radii.end());
    double max_binder_thickness = mat->binder_thickness_fraction * max_radii;
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

    simulator.add_periodic_boundary_condition('x', -box_side/2, box_side/2);
    simulator.add_periodic_boundary_condition('y', -box_side/2, box_side/2);

    std::cout << "Surfaces and periodic BCs generated"<<"\n";

    //************************************************
    // N.B. DOES NOT ALLOW FOR PARTICLES WITH BINDER !!
    // Generates the position of the particles with the function random_fill_box
    auto particle_positions = random_fill_box_periodic(-box_side/2, box_side/2, -box_side/2, box_side/2,
                                              -box_side/2, box_side/2, particle_radii, 0, "xy");
    std::cout << "Particle positions generated"<<"\n";
    //************************************************

    //Creates particles at particle positions
    for (std::size_t i = 0; i != particle_positions.size(); ++i)
    {
        simulator.create_particle(particle_radii[i], particle_positions[i], Vec3(0,0,0), mat);
    }

    std::chrono::duration<double> output_time { simulation_time/output_number};
    auto output1 = simulator.create_output(output_directory , output_time, "output");
    output1->print_particles = true;
    output1->print_fractured_particles = true;
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
    simulator.setup();
    EngineType::RunForTime run_for_time(simulator, simulation_time);

    bool restart = parameters.get_parameter<bool>("restart");


    for (unsigned level = 0; level != density_levels.size(); ++level)
    {
        double h_target = particle_volume/(box_side * box_side * density_levels[level]) - box_side/2;
        double surface_velocity = (h_target - top_surface->get_points()[0].z())/(simulation_time.count());
        std::cout << "Top surface position is : " << top_surface->get_points()[0].z() << "\n";
        std::cout << "Plate target is: " << h_target << "\n";
        std::cout << "Surface velocity is: " << surface_velocity << "\n";
        //===============Right Direction of surface disp===============
        top_surface->set_velocity(Vec3(0, 0, surface_velocity));
        //=============================================================

        run_for_time.reset(simulation_time);
        simulator.run(run_for_time);

//        else
//        {
        top_surface->set_velocity(Vec3(0,0,0));
        run_for_time.reset(simulation_time);
        simulator.run(run_for_time);
//        }
    }
    if (restart)
    {

        top_surface->set_velocity(Vec3(0,0,0));
        std::stringstream restart_file_name;
        restart_file_name << output_directory << "restart.res";
        simulator.write_restart_file(restart_file_name.str());
        std::cout << "Restarting simulation for unloading" << std::endl;
        auto unloading_simulator = Engine<ForceModel, ParticleType>(restart_file_name.str());
        std::cout << "New engine built" << std::endl;

        EngineType::RunForTime time_to_unload(unloading_simulator, simulation_time);
        std::cout << "before running unload sim" << std::endl;
        unloading_simulator.run(time_to_unload);

        time_to_unload.reset(simulation_time);
        unloading_simulator.run(time_to_unload);
    }
}