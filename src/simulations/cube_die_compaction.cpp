//
// Created by Axel on 2021-09-13.
//

#include <string>

#include "../engine/engine.h"
#include "../materials/elastic_ideal_plastic_material.h"
#include "../surfaces/point_surface.h"
#include "../utilities/file_reading_functions.h"
#include "../utilities/filling_functions.h"
#include "simulations.h"
#include "../contact_models/storakers_mesarovic_johnson.h"

void DEM::cube_die_compaction(const std::string& settings_file_name) {
    using namespace DEM;
    using ForceModel = StorakersMesarovicJohnson;
    using ParticleType = SphericalParticle<ForceModel>;

    // Defines the engine as an engine with StorakersMesarovicJohnson as force model and SphericalParticle as
    // particle type
    using EngineType = Engine<ForceModel, ParticleType>;

    // Makes it possible to write 1s for one second
    using namespace std::chrono_literals;

    // Reads the file defining different parameters which is passed as the second argument to the program
    // The typ of simulation is the first argument
    SimulationParameters parameters(settings_file_name);

    // Reads the parameter N from the simulation file which will be the number of particles
    auto N = parameters.get_parameter<std::size_t>("N");
    auto output_directory = parameters.get_parameter<std::string>("output_dir");
    auto particle_file = parameters.get_parameter<std::string>("radius_file");

    // Creates a DEM engine with a time step of one microsecond
    EngineType simulator(1us); //might have to be increased or decreased--> handpåläggning

    // Reads a vector with different relative densities where unloading will be made
    auto density_levels = parameters.get_vector<double>( "density_levels" );

    // Creates the material, Elastic - ideal plastic where the elastic part upon loading is assumed to be negligible
    // in the Storakers-Mesarovic-Johnson model
    auto material = simulator.create_material<ElasticIdealPlasticMaterial>(2630.); //assumed to be density of Alu in kg/m^3
    material->sY = parameters.get_parameter<double>("sY");
    material->E = parameters.get_parameter<double>("E");
    material->nu = parameters.get_parameter<double>("nu");

    material->mu = parameters.get_parameter<double>("mu");
    material->mu_wall = parameters.get_parameter<double>("mu_wall");
    material->kT = parameters.get_parameter<double>("kT");

    // Imports the initial relative density of the particles at filling
    auto particle_density_at_filling = parameters.get_parameter<double>("filling_density");


    // Importing the compaction time, the unloading velocity and the unloading time
    std::chrono::duration<double> compaction_time {parameters.get_parameter<double>("compaction_time")};
    auto unloading_velocity = parameters.get_parameter<double>("unloading_velocity");
    std::chrono::duration<double> unloading_time {parameters.get_parameter<double>("unloading_time")};



    // Read particle radii from file and creates a vector with N particles

    auto particle_radii = read_vector_from_file<double>(particle_file);
    particle_radii.assign(N,*particle_radii.begin());
//============================================================================================================


    // ================================================================================================================
    // Creating the particles
    // ================================================================================================================

    //  Calculating the volume of all particles
    double particle_volume = 0.;
    for(auto& r: particle_radii) {
        particle_volume += 4.*pi*r*r*r/3.;
    }
    std::cout << "Volume of simulated particles is " << particle_volume << "\n";

    //Creates the cube so that the volume corresponds to an initial density of filling_density

    auto box_side = pow(particle_volume/particle_density_at_filling, 1./3);
    std::cout << "box_side " <<box_side << "\n";

    //================
    //Creating the box
    //================

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


    std::cout << "Generating surfaces"<<"\n";
    // Creates all surfaces as point surfaces and writes out their normal-vector in the terminal
    auto top_surface = simulator.create_point_surface(top_points, true, "top_surface", false);
    std::cout << "Normal of top surface: " << top_surface->get_normal() << "\n";
    auto bottom_surface = simulator.create_point_surface(bottom_points, true, "bottom_surface", false);
    std::cout << "Normal of bottom surface: " << bottom_surface->get_normal() << "\n";
    auto front_surface = simulator.create_point_surface(front_points, true, "front_surface", false);
    std::cout << "Normal of front surface: " << front_surface->get_normal() << "\n";
    auto back_surface = simulator.create_point_surface(back_points, true, "back_surface", false);
    std::cout << "Normal of back surface: " << back_surface->get_normal() << "\n";
    auto right_surface = simulator.create_point_surface(right_points, true, "right_surface", false);
    std::cout << "Normal of right surface: " << right_surface->get_normal() << "\n";
    auto left_surface = simulator.create_point_surface(left_points, true, "left_surface", false);
    std::cout << "Normal of left surface: " << left_surface->get_normal() << "\n";

    std::cout << "Surfaces generated"<<"\n";

    //puts the rotation of the particles to 0
    simulator.set_rotation(false);
    std::cout << "Rotation locked"<<"\n";

    // Generates the position of the particles with the function random_fill_box
    auto particle_positions = random_fill_box(-box_side/2, box_side/2, -box_side/2, box_side/2,
                                              -box_side/2, box_side/2, particle_radii); //, material->bt); what material parameter is this?
    std::cout << "Particle positions generated"<<"\n";

    //Creates particles at particle positions
    for (std::size_t i = 0; i != particle_positions.size(); ++i) {
        simulator.create_particle(particle_radii[i], particle_positions[i], Vec3(0,0,0), material);
    }
    //creates an output for the filling, output every 5 millisecond
    auto output1 = simulator.create_output(output_directory , 0.005s, "output"); //time might be changed
    output1->print_particles = true;
    output1->print_kinetic_energy = true;
    output1->print_surface_positions = true;
    output1->print_surface_forces = true;


    //mass scaling
    simulator.set_mass_scale_factor(10.0); //how should this be chosen?

    //========================
    //Running the Simulation
    //========================
    simulator.setup();
    EngineType::RunForTime run_for_time(simulator, 0.1s);

    // calculating the distance required for prescribed density levels and giving each surface a velocity to reach a certain density
    for (unsigned level = 0; level != density_levels.size(); ++level) {
        double h_target = pow(particle_volume/density_levels[level], 1./3)/2;
        double surface_velocity = (h_target - top_surface->get_points()[0].z())/(compaction_time.count());
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


        //===============Wrong Direction of surface disp===============
        //top_surface->set_velocity(Vec3(0, 0, -surface_velocity));
        //bottom_surface->set_velocity(Vec3(0, 0, surface_velocity));
        //front_surface->set_velocity(Vec3(0, surface_velocity, 0));
        //back_surface->set_velocity(Vec3(0, -surface_velocity, 0));
        //left_surface->set_velocity(Vec3(surface_velocity, 0, 0));
        //right_surface->set_velocity(Vec3(-surface_velocity, 0, 0));
        //==============================================================
        run_for_time.reset(compaction_time);
        simulator.run(run_for_time);


        std::stringstream restart_file_name;
        restart_file_name << output_directory << "/restart_D=" << density_levels[level] << ".res";
        simulator.write_restart_file(restart_file_name.str());

        auto unloading_simulator = Engine<ForceModel, ParticleType>(restart_file_name.str());
        // Unloading the compact, setting an unloading velocity to all surfaces
        auto unload_top_surface = unloading_simulator.get_surface<EngineType::PointSurfacePointer>("top_surface");
        unload_top_surface->set_velocity(Vec3(0, 0, unloading_velocity));
        auto unload_bottom_surface = unloading_simulator.get_surface<EngineType::PointSurfacePointer>("bottom_surface");
        unload_bottom_surface->set_velocity(Vec3(0, 0, -unloading_velocity));
        auto unload_front_surface = unloading_simulator.get_surface<EngineType::PointSurfacePointer>("front_surface");
        unload_front_surface->set_velocity(Vec3(0, -unloading_velocity, 0));
        auto unload_back_surface = unloading_simulator.get_surface<EngineType::PointSurfacePointer>("back_surface");
        unload_back_surface->set_velocity(Vec3(0, unloading_velocity, 0));
        auto unload_left_surface = unloading_simulator.get_surface<EngineType::PointSurfacePointer>("left_surface");
        unload_left_surface->set_velocity(Vec3(-unloading_velocity, 0, 0));
        auto unload_right_surface = unloading_simulator.get_surface<EngineType::PointSurfacePointer>("right_surface");
        unload_right_surface->set_velocity(Vec3(unloading_velocity, 0, 0));

        auto compaction_output = unloading_simulator.get_output("output");
        unloading_simulator.remove_output(compaction_output);
        std::stringstream unloading_output_name;
        unloading_output_name << output_directory << "/unload_D=" << density_levels[level];
        auto output2 = unloading_simulator.create_output(unloading_output_name.str(), 0.001s);
        output2->print_particles = true;
        output2->print_kinetic_energy = true;
        output2->print_surface_positions = true;
        output2->print_surface_forces = true;
        restart_file_name << 1;
        unloading_simulator.write_restart_file(restart_file_name.str());
        EngineType::RunForTime time_to_unload(unloading_simulator, unloading_time);
        unloading_simulator.run(time_to_unload);
    }
}