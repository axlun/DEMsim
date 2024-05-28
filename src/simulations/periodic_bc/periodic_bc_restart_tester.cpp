//
// Created by Axel on 2024-03-18.
//

#include "../../engine/engine.h"
#include "../../materials/electrode_material.h"
#include "../../surfaces/point_surface.h"
#include "../../utilities/file_reading_functions.h"
#include "../../utilities/filling_functions.h"
#include "../simulations.h"
#include "../../contact_models/Positive_electrode/elastic_plastic_binder_elastic_plastic_particle.h"

void DEM::periodic_bc_restart_tester(const std::string& settings_file_name)
{
    using namespace DEM;
    using ForceModel = elastic_plastic_binder_elastic_plastic_particle;
    using ParticleType = SphericalParticle<ForceModel>;
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

    auto N = parameters.get_parameter<std::size_t>("N");
    std::cout << "N = " << N << std::endl;
    auto particle_file = parameters.get_parameter<std::string>("radius_file");

    // Reads a vector with different relative densities where unloading will be made
    auto density_levels = parameters.get_vector<double>("density_levels");

//    auto swell_rate = parameters.get_parameter<double>("swell_rate");
//    std::cout << "swell_rate = " << swell_rate << std::endl;
//
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
    // Imports the initial relative density of the particles at filling
    auto particle_density_at_filling = parameters.get_parameter<double>("filling_density");

    std::chrono::duration<double> simulation_time{parameters.get_parameter<double>("simulation_time")};


    // ================================================================================================================
    // Creating the particles
    // ================================================================================================================
    auto particle_radii = read_vector_from_file<double>(particle_file);
    particle_radii.assign(N, 1.0); // *****one particle size with radius 1*****
//    particle_radii.assign(N, *particle_radii.begin());
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
    for (const auto &r: particle_radii)
    {
        particle_volume += 4. / 3. * 3.1415 * r * r * r;
    }
    std::cout << "Volume of particles are: " << particle_volume << std::endl;
    //============================================================================================================

    //Creates the cube so that the volume corresponds to an initial density of filling_density

    auto box_side = pow(particle_volume / particle_density_at_filling, 1. / 3);
    double aspect_ratio = 1.5;
    auto box_height = box_side * aspect_ratio;
    box_side *= pow(1. / aspect_ratio, 0.5);
    std::cout << "box_height" << box_height<< "\n";
    std::cout << "box_side " << box_side << "\n";

    // ================================================================================================================
    //Creating the box
    // ================================================================================================================

    //Creates all the points that define the initial box
    auto p1 = Vec3(-box_side / 2, -box_side / 2, 0);
    auto p2 = Vec3(box_side / 2, -box_side / 2, 0);
    auto p3 = Vec3(box_side / 2, box_side / 2, 0);
    auto p4 = Vec3(-box_side / 2, box_side / 2, 0);
    auto p5 = Vec3(-box_side / 2, -box_side / 2, box_height);
    auto p6 = Vec3(box_side / 2, -box_side / 2, box_height);
    auto p7 = Vec3(box_side / 2, box_side / 2, box_height);
    auto p8 = Vec3(-box_side / 2, box_side / 2, box_height);

    //Saves the points corresponding to a surface in a vector
    std::vector <Vec3> front_points{p5, p6, p2, p1};
    std::vector <Vec3> back_points{p7, p8, p4, p3};
    std::vector <Vec3> right_points{p6, p7, p3, p2};
    std::vector <Vec3> left_points{p8, p5, p1, p4};
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

    simulator.add_periodic_boundary_condition('x', -box_side / 2, box_side / 2);
    simulator.add_periodic_boundary_condition('y', -box_side / 2, box_side / 2);

    std::cout << "Surfaces and periodic BCs generated" << "\n";

    //************************************************
    // N.B. DOES NOT ALLOW FOR PARTICLES WITH BINDER !!
    // Generates the position of the particles with the function random_fill_box
    auto particle_positions =
            random_fill_box_periodic(-box_side / 2, box_side / 2, -box_side / 2, box_side / 2,
                                     0, box_height, particle_radii, 0,"xy"); //, material->bt); what material parameter is this?
//    auto particle_positions = random_fill_box(-box_side/2, box_side/2, -box_side/2, box_side/2,
//                                              -box_side/2, box_side/2, particle_radii); //, material->bt); what material parameter is this?
    std::cout << "Particle positions generated" << "\n";
    //************************************************

    //Creates particles at particle positions
    for (std::size_t i = 0; i != particle_positions.size(); ++i)
    {
        simulator.create_particle(particle_radii[i], particle_positions[i], Vec3(0, 0, 0), mat);
    }
//    for (auto& p:simulator.get_particles())
//    {
//        p->set_swell_rate(swell_rate);
//    }

    std::chrono::duration<double> output_time{simulation_time / output_number};
    std::cout << "output_time = " << output_time.count() << std::endl;
    auto output1 = simulator.create_output(output_directory, output_time, "output");
    output1->print_particles = true;
    output1->print_kinetic_energy = true;
    output1->print_surface_positions = true;
    output1->print_surface_forces = true;
    output1->print_fabric_force_tensor = true;
    output1->print_contacts = true;
    output1->print_periodic_bc = true;
    output1->print_mirror_particles = true;

    bool restart = parameters.get_parameter<bool>("restart");
    simulator.set_gravity(Vec3(0, 0, -10));
    //========================
    //Running the Simulation
    //========================
    simulator.setup();
    EngineType::RunForTime run_for_time(simulator, simulation_time);

    if (restart)
    {
//        simulator.make_output();
        run_for_time.reset(simulation_time);
        simulator.run(run_for_time);
//        simulator.make_output();

        std::stringstream restart_file_name;
        restart_file_name << output_directory << "restart.res";
        simulator.write_restart_file(restart_file_name.str());
        std::cout << "Restarting simulation" << std::endl;
        auto restart_simulator = Engine<ForceModel, ParticleType>(restart_file_name.str());
        std::cout << "New engine built" << std::endl;
//        restart_simulator.make_output();
        EngineType::RunForTime time_after_restart(restart_simulator, simulation_time);
        restart_simulator.run(time_after_restart);
//        restart_simulator.make_output();
        std::stringstream restart_2_filename;
        restart_2_filename << output_directory << "restart_2.res";
        restart_simulator.write_restart_file(restart_2_filename.str());
    }
    else
    {
//        run_for_time.reset(simulation_time*0);
//        simulator.run(run_for_time);
//        simulator.make_output();

        run_for_time.reset(2*simulation_time);
        simulator.run(run_for_time);
//        simulator.make_output();
    }
}
//    for (unsigned level = 0; level != density_levels.size(); ++level)
//    {
//       if (restart)
//       {
//
//       }
//       else
//       {
//
//       }
//
//
//        double h_target = particle_volume/(box_side * box_side * density_levels[level]) - box_side/2;
//        double surface_velocity = (h_target - top_surface->get_points()[0].z())/(simulation_time.count());
//        std::cout << "Top surface position is : " << top_surface->get_points()[0].z() << "\n";
//        std::cout << "Plate target is: " << h_target << "\n";
//        std::cout << "Surface velocity is: " << surface_velocity << "\n";
//        ===============Right Direction of surface disp===============
//        top_surface->set_velocity(Vec3(0, 0, surface_velocity));
//        =============================================================
//
//        simulator.make_output();
//        run_for_time.reset(simulation_time);
//        simulator.run(run_for_time);
//        simulator.make_output();
//
//        top_surface->set_velocity(Vec3(0,0,0));
//        run_for_time.reset(simulation_time);
//        simulator.run(run_for_time);
//        if (restart)
//        {
//            std::stringstream restart_file_name;
//            restart_file_name << output_directory << "restart_D=" << density_levels[level] << ".res";
//            simulator.write_restart_file(restart_file_name.str());
//            std::cout << "Restarting simulation for unloading" << std::endl;
//            auto unloading_simulator = Engine<ForceModel, ParticleType>(restart_file_name.str());
//            std::cout << "New engine built" << std::endl;
//
//            std::stringstream restart_2_filename;
//            restart_2_filename << output_directory << "restart_2.res";
//            unloading_simulator.write_restart_file(restart_2_filename.str());
//        simulator.make_output();
//
//        auto compaction_output = unloading_simulator.get_output("output");
//        unloading_simulator.remove_output(compaction_output);
//            /*
//            std::stringstream unloading_output_name;
//            unloading_output_name << output_directory << "/unload_D=" << density_levels[level];
//            unloading_output_name << output_directory << "restart";
//            auto output2 = unloading_simulator.create_output(unloading_output_name.str(), output_time);
//            output2->print_particles = true;
//            output2->print_kinetic_energy = true;
//            output2->print_surface_positions = true;
//            output2->print_surface_forces = true;
//            output2->print_fabric_force_tensor = true;
//            output2->print_contacts = true;
//            output2->print_periodic_bc = true;
//            output2->print_mirror_particles = true;
//            unloading_simulator.make_output();
//            restart_file_name << 1;
//            unloading_simulator.write_restart_file(restart_file_name.str());
//            */
//
//             Unloading the compact, setting an unloading velocity to all surfaces
//            auto unload_top_surface =
//                    unloading_simulator.get_surface<EngineType::PointSurfacePointer>("top_surface");
//            unload_top_surface->set_velocity(Vec3(0, 0, 0));// -surface_velocity));
//            EngineType::RunForTime time_to_unload(unloading_simulator, simulation_time);
//            std::cout << "before running unload sim" << std::endl;
//            unloading_simulator.run(time_to_unload);
//        }
//        else
//        {
//            top_surface->set_velocity(Vec3(0,0,0));
//            run_for_time.reset(simulation_time);
//            simulator.run(run_for_time);
//        }
//    }
//}