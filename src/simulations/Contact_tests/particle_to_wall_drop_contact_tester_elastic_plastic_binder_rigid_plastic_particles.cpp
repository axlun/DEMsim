//
// Created by Axel on 2023-03-07.
//

#include <chrono>

#include "../simulations.h"

#include "../../engine/engine.h"
#include "../../particles/spherical_particle.h"
#include "../../surfaces/deformable_point_surface.h"
#include "../../engine/contact.h"
#include "../../utilities/file_reading_functions.h"
#include "../../contact_models/Positive_electrode/elastic_plastic_binder_rigid_plastic_particle.h"
#include "../../materials/electrode_material.h"
#include "../../materials/porous_electrode_material.h"

void DEM::particle_to_wall_drop_contact_tester_elastic_plastic_binder_rigid_plastic_particles(const std::string &settings_file_name)
{
    using namespace DEM;
    using ForceModel = elastic_plastic_binder_rigid_plastic_particle;
    using ParticleType = SphericalParticle<ForceModel>;
    using EngineType = Engine<ForceModel, ParticleType>;
    using SurfaceType = DeformablePointSurface<ForceModel,ParticleType>;
    using namespace std::chrono_literals;
    namespace fs = std::filesystem;


    SimulationParameters parameters{settings_file_name};
    auto output_directory = parameters.get_parameter<std::string>("output_dir");

    float time_step = parameters.get_parameter<double>("time_step")*1e-6;
    std::chrono::duration<double> time_step_us {time_step};
    std::cout << "Time step is:" << time_step_us.count()*1E6 << "µs\n";
    EngineType simulator(time_step_us); //orig  1E0

    auto mass_scaling = parameters.get_parameter<double>("mass_scaling");
    simulator.set_mass_scale_factor(mass_scaling); //Orig 1E2
    std::cout << "Mass scaling is:" << mass_scaling << "\n";

    auto radius = parameters.get_parameter<double>("R");
    auto height_factor = parameters.get_parameter<double>("height_factor");
    auto gravity = parameters.get_parameter<double>("gravity");
    auto Number_of_falls = parameters.get_parameter<double>("Number_of_falls");
    auto result_points = parameters.get_parameter<double>("result_points");

    auto mat = simulator.create_material<ElectrodeMaterial>(4800);
    mat->E = parameters.get_parameter<double>("E");
    mat->nu = parameters.get_parameter<double>("nu");
    mat->Ep=parameters.get_parameter<double>("Ep");
    mat->nup=parameters.get_parameter<double>("nup");
    mat->binder_radius_fraction=parameters.get_parameter<double>("binder_radius_fraction");
    mat->binder_thickness_fraction=parameters.get_parameter<double>("binder_thickness_fraction");
    mat->binder_stiffness_coefficient=parameters.get_parameter<double>("binder_stiffness_coefficient");
    mat->yield_displacement_coeff=parameters.get_parameter<double>("yield_displacement_coeff");
    mat->binder_yield_stress_=parameters.get_parameter<double>("binder_yield_stress_");
    mat->particle_yield_stress_=parameters.get_parameter<double>("particle_yield_stress_");
    mat->fraction_binder_contacts =parameters.get_parameter<double>("fraction_binder_contacts");
    mat->alpha_i = parameters.get_vector<double>("alpha_i");
    mat->tau_i =parameters.get_vector<double>("tau_i");
    mat->adhesive = false;
//    auto p1 = simulator.create_particle(radius,Vec3{0,0 , 0},Vec3{0,0,0}, mat);
//    auto p2 = simulator.create_particle(radius,Vec3{2*radius+(mat->binder_thickness_fraction*radius)*1.001 ,0 , 0},Vec3{0,0,0}, mat);
    double fall_height_center = 0.0;
    if (mat->fraction_binder_contacts > 0.0)
    {
        fall_height_center = height_factor * (radius + (mat->binder_thickness_fraction * radius) * 1.0001);
    }
    else {
        fall_height_center = height_factor * radius  * 1.0001 ;
    }
    std::cout << "Fall height: " << fall_height_center << "m\n";

    auto p1 = simulator.create_particle(radius,Vec3{0,0 , fall_height_center},Vec3{0,0,0}, mat);

    // Points for wall
    auto Wp1 = Vec3(2*radius, 2*radius, 0);
    auto Wp2 = Vec3(2*radius, -2*radius, 0);
    auto Wp3 = Vec3(-2*radius,  -2*radius, 0);
    auto Wp4 = Vec3(-2*radius,  2*radius, 0);
    std::vector<Vec3> wall_points{Wp1, Wp2, Wp3, Wp4};


    auto contact_surface = simulator.create_point_surface(wall_points, true,"contact_surface",true);
    std::cout << "Normal of contact surface: " << contact_surface->get_normal() << "\n";



//    simulator.set_mass_scale_factor(1E0);
    double bounding_box_stretch = 100 * radius*mat->binder_thickness_fraction;
    simulator.setup(bounding_box_stretch);
    simulator.set_gravity(Vec3(0, 0, -gravity));

    std::chrono::duration<double> fall_time(pow(2*(fall_height_center-radius)/gravity,.5));
    std::cout << "Fall time is: " << fall_time.count() <<"\n";
    std::chrono::duration<double> simulation_time = Number_of_falls * 2.0 * fall_time;
    std::cout << "Simulation time is: " << simulation_time.count() <<"\n";
    EngineType ::RunForTime run_for_time(simulator,simulation_time);


    std::chrono::duration<double> output_interval = simulation_time/result_points;
    std::cout << "Output interval is: " << output_interval.count() <<"\n";
    auto contact_output = simulator.create_output(output_directory,output_interval);
    contact_output->print_particles = true;
    contact_output->print_contacts = true;
    contact_output->print_fabric_force_tensor = true;
    contact_output->print_kinetic_energy = true;
    contact_output->print_surface_forces = true;
    contact_output->print_surface_positions = true;


    simulator.run(run_for_time);

/*    p1->set_velocity(Vec3{0,particle_velocity,0});
    p2->set_velocity(Vec3{0,-particle_velocity,0});
    auto particle_tangential_displacement =  radius*mat->binder_thickness_fraction*0.01; //Move particles 1% of the binder distance
    std::chrono::duration<double>particle_tangential_displacement_time(particle_tangential_displacement/particle_velocity);
    run_for_time.reset(particle_tangential_displacement_time);
    std::cout << "Tangential movement of particles"<< "\n";
    simulator.run(run_for_time);
*/
}