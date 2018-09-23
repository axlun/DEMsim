//
// Created by erolsson on 2018-09-18.
//

#include "simulations.h"

#include "engine.h"
#include "file_reading_functions.h"
#include "filling_functions.h"
#include "linear_contact_material.h"
#include "linear_stick_slip_model.h"
#include "vec3.h"


void DEM::gyratory_compaction(const std::string& settings_file_name){
    using namespace DEM;
    using ForceModel = LinearStickSlipModel;
    using ParticleType = SphericalParticle<ForceModel>;
    using EngineType = Engine<ForceModel, ParticleType>;

    using namespace std::chrono_literals;

    SimulationParameters parameters(settings_file_name);

    auto N = parameters.get<std::size_t>("N");
    auto output_directory = parameters.get<std::string>("output_dir");
    auto particle_file = parameters.get<std::string>("radius_file");
    auto filling_density = parameters.get<double>("filling_density");

    EngineType simulator(1us);

    auto material = simulator.create_material<LinearContactMaterial>(10e-9);
    material->k = 1000;
    material->mu=0.3;
    material->kT = 1000;
    material->mu_wall=0.3;

    // Read particle radii from file
    auto particle_radii = read_vector_from_file<double>(particle_file);
    particle_radii.assign(particle_radii.begin(), particle_radii.begin()+N);
    auto aspect_ratio_at_dense = parameters.get<double>("aspect_ratio_at_dense");
    double particle_volume = 0.;
    for(auto& r: particle_radii) {
        particle_volume += 4.*pi*r*r*r/3.;
    }

    std::cout << "Volume of simulated particles is " << particle_volume << "\n";
    auto cylinder_radius = pow(4*particle_volume/pi/aspect_ratio_at_dense, 1./3)/2;
    auto cylinder_height = 2*cylinder_radius*aspect_ratio_at_dense/filling_density;
    simulator.create_cylinder(cylinder_radius, Vec3(0, 0, 1), Vec3(1, 0, 0), 2, true, true);
    std::cout << "The simulated cylinder has a radius of " << cylinder_radius << " and a height of "
              << cylinder_height << "\n";
    auto particle_positions = random_fill_cylinder(0, cylinder_height, cylinder_radius, particle_radii);

    for (std::size_t i=0; i != particle_positions.size(); ++i) {
        simulator.create_particle(particle_radii[i], particle_positions[i], Vec3(0,0,0), material);
    }


    // Creating The bottom plate surface
    Vec3 p1(-cylinder_radius, -cylinder_radius, 0.);
    Vec3 p2(-cylinder_radius,  cylinder_radius, 0.);
    Vec3 p3(cylinder_radius,   cylinder_radius, 0.);
    Vec3 p4(cylinder_radius,  -cylinder_radius, 0.);
    std::vector<Vec3> bottom_points{p4, p3, p2, p1};
    auto bottom_surface = simulator.create_point_surface(bottom_points, true);
    std::cout << "Normal of bottom surface is " << bottom_surface->get_normal() << std::endl;

    // Creating The top plate surface
    Vec3 p5(-cylinder_radius, -cylinder_radius, cylinder_height);
    Vec3 p6(-cylinder_radius,  cylinder_radius, cylinder_height);
    Vec3 p7(cylinder_radius,   cylinder_radius, cylinder_height);
    Vec3 p8(cylinder_radius,  -cylinder_radius, cylinder_height);
    std::vector<Vec3> top_points{p5, p6, p7, p8};
    auto top_surface = simulator.create_point_surface(top_points, true);
    std::cout << "Normal of top surface is " << top_surface->get_normal() << std::endl;

    auto output1 = simulator.create_output(output_directory, 0.01s);
    output1->print_particles = true;
    output1->print_kinetic_energy = true;

    simulator.set_gravity(Vec3(0, 0, -9820));
    simulator.setup();
    EngineType::RunForTime run_for_time(simulator, 0.05s);
    simulator.run(run_for_time);
}