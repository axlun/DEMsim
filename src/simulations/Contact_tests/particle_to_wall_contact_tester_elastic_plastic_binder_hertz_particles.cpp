//
// Created by Axel on 2022-08-30.
//

#include <chrono>

#include "../simulations.h"

#include "../../engine/engine.h"
#include "../../particles/spherical_particle.h"
#include "../../surfaces/deformable_point_surface.h"
#include "../../engine/contact.h"
#include "../../utilities/file_reading_functions.h"
#include "../../contact_models/viscoelastic.h"
#include "../../contact_models/porous_electrode_contact.h"
#include "../../contact_models/viscoelastic_binder_with_deformable_particles.h"
#include "../../contact_models/elastic_perfect_plastic.h"
//#include "../../contact_models/Positive_electrode/elastic_plastic_binder_rigid_perfect_plastic_particle_OLD.h"
#include "../../contact_models/Positive_electrode/elastic_plastic_binder_hertz_plastic_particle.h"
#include "../../materials/electrode_material.h"
#include "../../materials/porous_electrode_material.h"

void DEM::particle_to_wall_contact_tester_elastic_plastic_binder_hertz_particles(const std::string &settings_file_name)
{
    using namespace DEM;
    using ForceModel = elastic_plastic_binder_hertz_plastic_particle;
    using ParticleType = SphericalParticle<ForceModel>;
    using EngineType = Engine<ForceModel, ParticleType>;
    using SurfaceType = DeformablePointSurface<ForceModel,ParticleType>;
    using namespace std::chrono_literals;
    namespace fs = std::filesystem;

    EngineType simulator(1E-0us); //orig  1E0
    SimulationParameters parameters{settings_file_name};

    auto mat = simulator.create_material<ElectrodeMaterial>(4800);
//    using Material = ElectrodeMaterial;
//    Material mat {0, 4800.};


    auto radius = parameters.get_parameter<double>("R");
    auto height_factor = parameters.get_parameter<double>("height_factor");
    auto gravity = parameters.get_parameter<double>("gravity");
    auto result_points = parameters.get_parameter<double>("result_points");
    auto output_directory = parameters.get_parameter<std::string>("output_dir");
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
    auto fall_height = height_factor*(radius+(mat->binder_thickness_fraction*radius)*1.0001);
    auto p1 = simulator.create_particle(radius,Vec3{0,0 , fall_height},Vec3{0,0,0}, mat);

    // Points for wall
    auto Wp1 = Vec3(2*radius, 2*radius, 0);
    auto Wp2 = Vec3(2*radius, -2*radius, 0);
    auto Wp3 = Vec3(-2*radius,  -2*radius, 0);
    auto Wp4 = Vec3(-2*radius,  2*radius, 0);
    std::vector<Vec3> wall_points{Wp1, Wp2, Wp3, Wp4};


    auto contact_surface = simulator.create_point_surface(wall_points, true,"contact_surface",true);
    std::cout << "Normal of contact surface: " << contact_surface->get_normal() << "\n";



    simulator.set_mass_scale_factor(1E0);
    double bounding_box_stretch = radius*mat->binder_thickness_fraction;
    simulator.setup(bounding_box_stretch);
    simulator.set_gravity(Vec3(0, 0, -gravity));

    std::chrono::duration<double> fall_time(pow(2*fall_height/gravity,.5)); //Divide by 2 as both particles are moving
    EngineType ::RunForTime run_for_time(simulator,fall_time);
    std::cout << "Fall time is: " << fall_time.count() <<"\n";


    std::chrono::duration<double> output_interval = fall_time/result_points;
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