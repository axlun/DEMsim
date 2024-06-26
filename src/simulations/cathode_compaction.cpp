//
// Created by Axel on 2021-09-27.
//

#include "simulations.h"

#include <vector>

#include "../engine/engine.h"
#include "../contact_models/viscoelastic.h"
#include "../materials/electrode_material.h"
#include "../utilities/file_reading_functions.h"
#include "../utilities/filling_functions.h"



void DEM::cathode_compaction(const std::string& settings_file_name) {
    using namespace DEM;
    using ForceModel = Viscoelastic;
    using ParticleType = SphericalParticle<ForceModel>;
    using EngineType = Engine<ForceModel, ParticleType>;
    using namespace std::chrono_literals;
    SimulationParameters parameters(settings_file_name);
    auto N = parameters.get_parameter<double>("N");
    auto output_directory = parameters.get_parameter<std::string>("output_dir");
    auto particle_file = parameters.get_parameter<std::string>("radius_file");

    EngineType simulator(1us);

    auto mat = simulator.create_material<ElectrodeMaterial>(4800);
    mat->E = parameters.get_parameter<double>("E");
    //mat->kT = parameters.get_parameter<double>("kT");
    mat->Ep = parameters.get_parameter<double>("Ep");
    mat->nu = parameters.get_parameter<double>("nu");
    mat->fraction_binder_contacts = parameters.get_parameter<double>("fraction_binder_contacts");
    mat->binder_radius_fraction = parameters.get_parameter<double>("binder_radius_fraction");
    mat->nup = parameters.get_parameter<double>("nup");
    mat->mu = parameters.get_parameter<double>("mu");
    mat->mu_wall = parameters.get_parameter<double>("mu_wall");
    mat->tau_i = parameters.get_vector<double>("tau_i");
    mat->alpha_i = parameters.get_vector<double>("alpha_i");
    mat->bt = parameters.get_parameter<double>("bt");
    mat->mu_binder = parameters.get_parameter<double>("mu_binder");
    auto particle_density_at_filling = parameters.get_parameter<double>("filling_density");
    auto particle_density_at_cube = parameters.get_parameter<double>("particle_density_at_cube");
    mat->active_particle_height = parameters.get_parameter<double>("active_particle_height");


    auto particle_radii = read_vector_from_file<double>(particle_file);
    particle_radii.assign(N, *particle_radii.begin());

    double particle_volume = 0.;
    for (const auto &r: particle_radii)
    {
        particle_volume += 4.*pi*r*r*r/3.;
    }
    std::cout << "Volume of simulated particles is " << particle_volume << "\n";
    auto box_side = pow(particle_volume/particle_density_at_cube, 1./3); //Boxens sida beroende av volymen av partiklarna
    // och relativa densiteten när partiklarna är kompakterade
    std::cout << "box_side " <<box_side << "\n";
    auto box_height = particle_density_at_cube*box_side/particle_density_at_filling;
    std::cout << "box_height " <<box_height << "\n";


    auto p1 = Vec3(-box_side/2, -box_side/2, 0);
    auto p2 = Vec3( box_side/2, -box_side/2, 0);
    auto p3 = Vec3( box_side/2,  box_side/2, 0);
    auto p4 = Vec3(-box_side/2,  box_side/2, 0);
    auto p5 = Vec3(-box_side/2, -box_side/2, box_height);
    auto p6 = Vec3( box_side/2, -box_side/2, box_height);
    auto p7 = Vec3( box_side/2,  box_side/2, box_height);
    auto p8 = Vec3(-box_side/2,  box_side/2, box_height);
    std::vector<Vec3> bottom_points{p1, p2, p3, p4};
    std::vector<Vec3> top_points{p8, p7, p6, p5};

    auto particle_positions = random_fill_box(-box_side/2, box_side/2, -box_side/2, box_side/2,
                                              0, box_height, particle_radii, mat->bt);

    auto deformable_surface = simulator.create_deformable_point_surface(bottom_points, true);
    auto top_surface = simulator.create_point_surface(top_points, true, "top_plate", false);

    std::cout << "Normal of top surface: " << top_surface->get_normal() << "\n";

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

    mat->adhesive = false;

    simulator.set_gravity(Vec3(0, 0, -9.82)); //comment
    simulator.set_mass_scale_factor(10.0);
    simulator.setup(1.01*mat->bt); //vad är materialparametern bt?
    simulator.set_rotation(false);
    EngineType::RunForTime run_for_time(simulator, 0.1s);
    simulator.run(run_for_time);
    EngineType::ParticleVelocityLess max_velocity (simulator, 0.1, 0.001s); //varför valdes denhär tiden?
    simulator.run(max_velocity);

    std::cout<<"beginning of compaction"<< std::endl;
    auto bbox = simulator.get_bounding_box();
    double h = bbox[5];
    top_surface->move(-Vec3(0, 0, box_height - h-1.01*mat->bt), Vec3(0, 0, 0));
    std::cout<<"h"<< h<< std::endl;
    double surface_velocity = 0.043; //varför denhär ythastigheten?
    mat-> adhesive = true;
    top_surface->set_velocity(Vec3(0, 0, 0.-surface_velocity));
    simulator.set_mass_scale_factor(10.0);
    std::chrono::duration<double> compaction_time {((h - mat->active_particle_height) / surface_velocity)};
    run_for_time.reset(compaction_time);
    simulator.run(run_for_time);

}