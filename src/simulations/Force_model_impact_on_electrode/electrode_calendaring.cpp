//
// Created by Axel on 2021-12-06.
//
#include "../simulations.h"

#include <vector>

#include "../../engine/engine.h"
#include "../../contact_models/viscoelastic.h"
#include "../../contact_models/Binder_behavour_investigation/viscoelastic_binder_El_Pl_particles.h"
#include "../../materials/electrode_material.h"
#include "../../utilities/file_reading_functions.h"
#include "../../utilities/filling_functions.h"

void DEM::electrode_calendaring(const std::string& settings_file_name) {
    using namespace DEM;
    using ForceModel = viscoelastic_binder_El_Pl_particles;
    using ParticleType = SphericalParticle<ForceModel>;
    using EngineType = Engine<ForceModel, ParticleType>;
    using namespace std::chrono_literals;

    SimulationParameters parameters(settings_file_name);
    auto output_directory = parameters.get_parameter<std::string>("output_dir");

    EngineType simulator(1us);

    auto N = parameters.get_parameter<double>("N");
    auto particle_file = parameters.get_parameter<std::string>("radius_file");
    auto mat = simulator.create_material<ElectrodeMaterial>(4800);
    mat->E = parameters.get_parameter<double>("E");
    mat->Ep = parameters.get_parameter<double>("Ep");
    mat->nup = parameters.get_parameter<double>("nup");
    mat->yield_displacement_coeff = parameters.get_parameter<double>("yield_displacement_coeff");
//    mat->mu = parameters.get_parameter<double>("mu");
//    mat->mu_wall = parameters.get_parameter<double>("mu_wall");
    mat->tau_i = parameters.get_vector<double>("tau_i");
    mat->alpha_i = parameters.get_vector<double>("alpha_i");
    mat->binder_thickness_fraction = parameters.get_parameter<double>("binder_thickness_fraction");
    mat->binder_radius_fraction = parameters.get_parameter<double>("binder_radius_fraction");
    mat->binder_stiffness_coefficient = parameters.get_parameter<double>("binder_stiffness_coefficient");
    mat->fraction_binder_contacts = parameters.get_parameter<double>("fraction_binder_contacts");
//  mat->mu_binder = parameters.get_parameter<double>("mu_binder");
    auto particle_density_at_filling = parameters.get_parameter<double>("filling_density");
    auto particle_density_at_cube = parameters.get_parameter<double>("particle_density_at_cube");
    mat->active_particle_height = parameters.get_parameter<double>("active_particle_height");

    auto particle_radii = read_vector_from_file<double>(particle_file);
    particle_radii.assign(N, *particle_radii.begin());

    double particle_volume = 0.;
    for (const auto &r: particle_radii) {
        particle_volume += 4. / 3. * pi * r * r * r;
    }
    std::cout << "Volume of particles is " << particle_volume << "\n";


    auto box_side = pow(particle_volume / particle_density_at_cube, 1. / 3.);
    std::cout << "box_side " << box_side << "\n";
    auto box_height = particle_density_at_cube * box_side / particle_density_at_filling;
    box_height = 1; //testing box height
    std::cout << "box_height " << box_height << "\n";

    auto p1 = Vec3(-box_side / 2, -box_side / 2, 0);
    auto p2 = Vec3(box_side / 2, -box_side / 2, 0);
    auto p3 = Vec3(box_side / 2, box_side / 2, 0);
    auto p4 = Vec3(-box_side / 2, box_side / 2, 0);
    auto p5 = Vec3(-box_side / 2, -box_side / 2, box_height);
    auto p6 = Vec3(box_side / 2, -box_side / 2, box_height);
    auto p7 = Vec3(box_side / 2, box_side / 2, box_height);
    auto p8 = Vec3(-box_side / 2, box_side / 2, box_height);
    std::vector <Vec3> bottom_points{p1, p2, p3, p4};
    std::vector <Vec3> top_points{p8, p7, p6, p5};

    //auto particle_positions = random_fill_box(-box_side / 2, box_side / 2, -box_side / 2, box_side / 2,
    //                                         0, box_height, particle_radii, mat->bt);
    DEM::Vec3 position = {0,0,.01};
    std::vector<Vec3> particle_positions = {};
    particle_positions.push_back(position);
    std::cout << "Particle positions: \n";
    for (const auto &i: particle_positions) {
        std::cout << i << "\n";
    }
    //auto deformable_surface = simulator.create_deformable_point_surface(bottom_points, true);
    auto deformable_surface = simulator.create_deformable_point_surface(bottom_points,"bottom_plate", true);
    auto top_surface = simulator.create_point_surface(top_points, true, "top_plate", false);

    std::cout << "Normal of top surface: " << top_surface->get_normal() << "\n";
    std::cout << "Normal of bottom surface: " << deformable_surface->get_normal() << "\n";

    for (std::size_t i = 0; i != particle_positions.size(); ++i) {
        simulator.create_particle(particle_radii[i], particle_positions[i], Vec3(0,0,0), mat);
    }
    auto filling_output = simulator.create_output(output_directory , 0.005s);
    filling_output->print_particles = true;
    filling_output->print_kinetic_energy = true;
    filling_output->print_surface_positions = true;
    filling_output->print_surface_forces = true;
    filling_output->print_contacts = true;
    filling_output->print_periodic_bc = true;
    filling_output->print_mirror_particles = true;
    filling_output->print_fabric_force_tensor=true;

    simulator.add_periodic_boundary_condition('x', -box_side/2, box_side/2);
    simulator.add_periodic_boundary_condition('y', -box_side/2, box_side/2);


    //Initial packing of particles, let particles fall with gravity
    mat->adhesive = true; //No adhesion of particles when initial packing
    simulator.set_gravity(Vec3(0, 0, -9.82)); //Use gravity for inital packing of particles
    simulator.set_mass_scale_factor(10.0);
    double max_radii = *max_element(particle_radii.begin(), particle_radii.end()); //find the biggest binder thickness to use when defineing a bounding box
    double max_binder_thickness = mat->binder_thickness_fraction*max_radii;
    //std::cout << "max_binder_thickness: "<< max_binder_thickness <<"\n";
    simulator.setup(1.01*max_binder_thickness); //Size of box for detecting contacts between particles

    simulator.set_rotation(false);
    std::cout << "Initialize natural particle packing \n";


    //Run for 0.1s and then run untill max_velocity of the paricles is 0.1 m/s and check it every 0.02s
    EngineType::RunForTime run_for_time(simulator, 0.1s);
    simulator.run(run_for_time);
    EngineType::ParticleVelocityLess max_velocity (simulator, 0.1, 0.02s); //Why this update time?
    simulator.run(max_velocity);

    //Calendaring process
    std::cout << "Initialize calendaring process \n";
    auto bbox = simulator.get_bounding_box(); //get the XYZ max/min that contain all particles
    double h = bbox[5]; //height of uppermost particle (Z-max)
    std::cout<<"h"<< h<< std::endl;
    top_surface->move(-Vec3(0, 0, box_height - h-1.01*max_binder_thickness), Vec3(0, 0, 0)); //Move topsurface to uppermost partile+binder thickness
    auto surface_velocity =parameters.get_parameter<double>("calendaring_surface_velocity"); //How was this chosen?
    mat->adhesive = true;
    top_surface->set_velocity(Vec3(0,0,0.-surface_velocity));
    std::chrono::duration<double> compaction_time {((h - mat->active_particle_height) / surface_velocity)};
    run_for_time.reset(compaction_time);
    simulator.run(run_for_time);
    simulator.write_restart_file(output_directory + "/compact_restart_file.res");

    std::cout<<"beginning of unloading"<< std::endl;
    top_surface->set_velocity(Vec3(0, 0, surface_velocity));
    EngineType::SurfaceNormalForceLess zero_force(top_surface, 0.);
    // simulator.set_rotation(false);
    simulator.run(zero_force);
}










