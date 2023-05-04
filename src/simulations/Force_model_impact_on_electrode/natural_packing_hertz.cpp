//
// Created by Axel on 2021-12-06.
//
#include "../simulations.h"

#include <vector>

#include "../../engine/engine.h"
#include "../../contact_models/viscoelastic.h"
#include "../../contact_models/Positive_electrode/elastic_plastic_binder_hertz_plastic_particle.h"
#include "../../materials/electrode_material.h"
#include "../../utilities/file_reading_functions.h"
#include "../../utilities/filling_functions.h"

void DEM::natural_packing_hertz(const std::string& settings_file_name) {
    using namespace DEM;
    using ForceModel = elastic_plastic_binder_hertz_plastic_particle;
    using ParticleType = SphericalParticle<ForceModel>;
    using EngineType = Engine<ForceModel, ParticleType>;
    using namespace std::chrono_literals;

    SimulationParameters parameters(settings_file_name);
    auto output_directory = parameters.get_parameter<std::string>("output_dir");
    float time_step = parameters.get_parameter<double>("time_step")*1e-6;
    std::chrono::duration<double> time_step_us {time_step};
//    std::chrono::duration<double> time_step_us {std::chrono::microseconds (1*time_step)};
    std::cout << "Time step is:" << time_step_us.count() << "\n";
    EngineType simulator(time_step_us); //orig  1E0

    auto N = parameters.get_parameter<double>("N"); //Number of particles
    auto particle_file = parameters.get_parameter<std::string>("radius_file");
    auto mat = simulator.create_material<ElectrodeMaterial>(4800);
    auto gravity = parameters.get_parameter<double>("gravity");
    auto mass_scaling = parameters.get_parameter<double>("mass_scaling");
    auto outputs = parameters.get_parameter<double>("outputs");
    mat->E = parameters.get_parameter<double>("E");
    mat->Ep = parameters.get_parameter<double>("Ep");
    mat->nu = parameters.get_parameter<double>("nu");
    mat->nup = parameters.get_parameter<double>("nup");
    mat->rhop = parameters.get_parameter<double>("rhop");
    mat->mu = parameters.get_parameter<double>("mu");
    mat->mu_wall = parameters.get_parameter<double>("mu_wall");
    auto rho_al = parameters.get_parameter<double>("rho_al"); //density of active layer
    auto mass_ratio_particles = parameters.get_parameter<double>("mass_ratio_particles"); //ratio between particle mass and total mass
    mat->tau_i = parameters.get_vector<double>("tau_i");
    mat->alpha_i = parameters.get_vector<double>("alpha_i");
//    mat->particle_yield_stress_ = parameters.get_parameter<double>("particle_yield_stress_");
    mat->binder_radius_fraction = parameters.get_parameter<double>("binder_radius_fraction");
    mat->binder_thickness_fraction = parameters.get_parameter<double>("binder_thickness_fraction");
    mat->binder_yield_stress_ = parameters.get_parameter<double>("binder_yield_stress_");
    mat->binder_stiffness_coefficient = parameters.get_parameter<double>("binder_stiffness_coefficient");
    mat->fraction_binder_contacts = parameters.get_parameter<double>("fraction_binder_contacts");
    auto particle_density_at_filling = parameters.get_parameter<double>("filling_density");
    mat->active_particle_height = parameters.get_parameter<double>("active_particle_height");
    mat->yield_displacement_coeff = parameters.get_parameter<double>("yield_displacement_coeff");

    auto particle_radii = read_vector_from_file<double>(particle_file);
    if (particle_radii.size() ==1) {
        particle_radii.assign(N, *particle_radii.begin());
    } else{
        particle_radii.assign(particle_radii.begin(), particle_radii.begin()+N);
        std::sort(particle_radii.rbegin(), particle_radii.rend());
    }
    double max_radii = *max_element(particle_radii.begin(), particle_radii.end()); //find the biggest binder thickness to use when defining a bounding box
    double max_binder_thickness = mat->binder_thickness_fraction*max_radii;

    double particle_volume = 0.;
    for (const auto &r: particle_radii) {
        particle_volume += 4. / 3. * pi * r * r * r;
        //std::cout << "Particle radii " << r << "\n";
    }
    std::cout << "Volume of particles is " << particle_volume << "\n";

    auto box_side = pow(particle_volume*mat->rhop/mat->active_particle_height/rho_al/mass_ratio_particles, 1./2.);
    std::cout << "box_side " << box_side << "\n";
    auto box_height = particle_volume/particle_density_at_filling/pow(box_side,2);
//    box_height = 1; //testing box height
    std::cout << "box_height " << box_height << "\n";

    auto p1 = Vec3(-box_side / 2, -box_side / 2, 0);
    auto p2 = Vec3(box_side / 2, -box_side / 2, 0);
    auto p3 = Vec3(box_side / 2, box_side / 2, 0);
    auto p4 = Vec3(-box_side / 2, box_side / 2, 0);
    auto p5 = Vec3(-box_side / 2, -box_side / 2, box_height);
    auto p6 = Vec3(box_side / 2, -box_side / 2, box_height);
    auto p7 = Vec3(box_side / 2, box_side / 2, box_height);
    auto p8 = Vec3(-box_side / 2, box_side / 2, box_height);

/*
    auto p1_stiff = Vec3(-.75*box_side / 2, -.75*box_side / 2, 0);
    auto p2_stiff = Vec3(.75*box_side / 2, -.75*box_side / 2, 0);
    auto p3_stiff = Vec3(.75*box_side / 2, .75*box_side / 2, 0);
    auto p4_stiff = Vec3(-.75*box_side / 2, .75*box_side / 2, 0);
    auto p5_stiff = Vec3(-.75*box_side / 2, -.75*box_side / 2, box_height);
    auto p6_stiff = Vec3(.75*box_side / 2, -.75*box_side / 2, box_height);
    auto p7_stiff = Vec3(.75*box_side / 2, .75*box_side / 2, box_height);
    auto p8_stiff = Vec3(-.75*box_side / 2, .75*box_side / 2, box_height);
*/
    std::vector <Vec3> bottom_points{p1, p2, p3, p4};
    std::vector <Vec3> top_points{p8, p7, p6, p5};

    //side panels for non periodic BCs
    std::vector <Vec3> side_1{p2, p3, p7, p6};
    std::vector <Vec3> side_2{p3, p4, p8, p7};
    std::vector <Vec3> side_3{p4, p1, p5, p8};
    std::vector <Vec3> side_4{p1, p2, p6, p5};

    auto particle_positions = random_fill_box(-box_side / 2, box_side / 2, -box_side / 2, box_side / 2,
                                             0, 0+box_height, particle_radii, max_binder_thickness);

// ******Specify position of two particles ******
//    DEM::Vec3 position1 = {0,0,3.9E-2};
//    DEM::Vec3 position2 = {0,0,40E-2};
//    std::vector<Vec3> particle_positions = {};
//    particle_positions.push_back(position1);
//    particle_positions.push_back(position2);
//    std::cout << "Particle positions: \n";
//    for (const auto &i: particle_positions) {
//        std::cout << i << "\n";
//    }

    auto deformable_surface = simulator.create_deformable_point_surface(bottom_points,"bottom_plate", true);
    auto top_surface = simulator.create_point_surface(top_points, true, "top_plate", false);

    //Side surface for non periodic BCs
    auto side1_surface = simulator.create_point_surface(side_1, true, "side1_plate", false);
    auto side2_surface = simulator.create_point_surface(side_2, true, "side2_plate", false);
    auto side3_surface = simulator.create_point_surface(side_3, true, "side3_plate", false);
    auto side4_surface = simulator.create_point_surface(side_4, true, "side4_plate", false);

    std::cout << "Normal of top surface: " << top_surface->get_normal() << "\n";
    std::cout << "Normal of bottom surface: " << deformable_surface->get_normal() << "\n";

    std::cout << "Normal of side surface 1: " << side1_surface->get_normal() << "\n";
    std::cout << "Normal of side surface 2: " << side2_surface->get_normal() << "\n";
    std::cout << "Normal of side surface 3: " << side3_surface->get_normal() << "\n";
    std::cout << "Normal of side surface 4: " << side4_surface->get_normal() << "\n";



    std::chrono::duration<double> packing_time {1 * pow((2*box_height/gravity),0.5)};
    std::cout<<"Gravity is: "<< gravity << '\n';
    std::cout<<"Box height is: "<< box_height << '\n';
    std::cout << "Natural packing time is: "<< packing_time.count() << '\n';
    std::cout << "Natural packing time is: "<< packing_time.count() << '\n';


    double packing_vel = packing_time.count()*gravity/10;
    std::cout << "packing_vel is: "<< packing_vel << '\n';

    std::chrono::duration<double> packing_time_constant_vel {(1.0/3.0)*pow(1*packing_time.count(),2)*gravity/2.0/packing_vel};
    std::cout << "Constnat velocity packing time is: "<< packing_time_constant_vel.count() << '\n';


    for (std::size_t i = 0; i != particle_positions.size(); ++i) {
        simulator.create_particle(particle_radii[i], particle_positions[i], Vec3(0,0,-packing_vel), mat);
    }

//=======================================PERIODIC BC:S===============================================================
//    simulator.add_periodic_boundary_condition('x', -box_side/2, box_side/2);
//    simulator.add_periodic_boundary_condition('y', -box_side/2, box_side/2);
//=====================================================================================================================
//new comment
    //Initial packing of particles, let particles fall with gravity
    mat->adhesive = false; //No adhesion of particles when initial packing
    simulator.set_mass_scale_factor(mass_scaling); //Orig 1E2
    std::cout << "max_binder_thickness: "<< max_binder_thickness <<"\n";
    simulator.setup(1.01*max_binder_thickness); //Size of box for detecting contacts between particles
    simulator.set_rotation(false);


// New packing method
    simulator.set_gravity(Vec3(0, 0, -0)); //Use gravity for initial packing of particles
    EngineType::RunForTime run_for_time(simulator, packing_time_constant_vel);


    std::chrono::duration<double> output_interval {packing_time_constant_vel.count()/outputs};
    std::cout << "Output interval: " << output_interval.count() << "\n";
    auto filling_output = simulator.create_output(output_directory , output_interval);
    filling_output->print_particles = true;
    filling_output->print_kinetic_energy = true;
    filling_output->print_surface_positions = true;
    filling_output->print_surface_forces = true;
    filling_output->print_contacts = true;
    filling_output->print_periodic_bc = true;
    filling_output->print_mirror_particles = true;
    filling_output->print_fabric_force_tensor=true;
    std::cout << "****************Initialize natural particle packing**************** \n";
    simulator.run(run_for_time);

    std::cout << "****************Turn on gravity**************** \n";

    simulator.set_gravity(Vec3(0,0,-gravity));
    run_for_time.reset(2*packing_time_constant_vel);
    simulator.run(run_for_time);

    std::cout << "****************Stop all particles and let come to rest**************** \n";
    for (auto& p: simulator.get_particles())
    {
        p->set_velocity(Vec3(0,0,0));
    }
    run_for_time.reset(2*packing_time_constant_vel);
    simulator.run(run_for_time);

    /*
    std::cout << "****************Initialize pre-calendering process**************** \n";
//  Move surface to uppermost particle
    auto bbox = simulator.get_bounding_box(); //get the XYZ max/min that contain all particles
    double h_1 = bbox[5]; //height of uppermost particle (Z-max)
    if (h_1<3.1*calendering_height){
        h_1=3.1*calendering_height;
        std::cout<<"Heigh of uppermoast particle lower then 2.2 h_al: "<< h_1<< std::endl;
    }
    else{
        std::cout<<"Heigh of uppermoast particle: "<< h_1<< std::endl;
    }
    top_surface->move(-Vec3(0, 0, box_height - h_1-1.01*max_binder_thickness), Vec3(0, 0, 0)); //Move top surface to uppermost partile+binder thickness

    top_surface->set_velocity(Vec3(0,0,0.-2*surface_velocity));
    std::chrono::duration<double> compaction_time_pre_cal {((h_1+1.01*max_binder_thickness - calendering_height*3) / (2*surface_velocity))};

    std::cout<<"Compaction time: "<< compaction_time_pre_cal.count()<< std::endl;
    run_for_time.reset(compaction_time_pre_cal);
    simulator.run(run_for_time);

    top_surface->set_velocity(Vec3(0,0,0));
    simulator.set_gravity(Vec3(0, 0, -1E1)); //Use gravity for initial packing of particles

    run_for_time.reset(4e0s);
    simulator.run(run_for_time);

//    EngineType::ParticleVelocityLess max_velocity_2 (simulator, 10, 0.04s);//2.5, 0.04s); // max_vel = 0.5
//    simulator.run(max_velocity_2);
// | ***********************************Move stiff surfaces and initiate periodic BCs***********************************
    std::cout << "****************Wall removal**************** \n";

    // Stop all the particles
    for (auto& p: simulator.get_particles())
    {
        p->set_velocity(Vec3(0,0,0));
    }

// |=====================MOVE THE STIFF SURFACE TO INITIATE THE PERIODIC BC:S ==========================================
    side1_surface->move(Vec3(5*box_side,0,0), Vec3(0, 0, 0));
    side2_surface->move(Vec3(0,5*box_side,0), Vec3(0,0,0));
    side3_surface->move(-Vec3(5*box_side,0,0), Vec3(0,0,0));
    side4_surface->move(-Vec3(0,5*box_side,0), Vec3(0,0,0));
//|=====================================================================================================================

    run_for_time.reset(2e0s);
    simulator.run(run_for_time);
    std::cout << "****************Adhesive on**************** \n";

    mat->adhesive = true; // Activate adhesion before calendering starts
    run_for_time.reset(2e0s);
    simulator.run(run_for_time);
//    max_velocity_2.set_new_value(3);//1.5);
//    simulator.run(max_velocity_2);
    std::cout << "****************Turn off gravity**************** \n";

    //turn off gravity
    simulator.set_gravity(Vec3(0, 0, 0));
    run_for_time.reset(1e0s);
    simulator.run(run_for_time);

    // Stop all the particles
    for (auto& p: simulator.get_particles())
    {
        p->set_velocity(Vec3(0,0,0));
    }
    run_for_time.reset(1e0s);
    simulator.run(run_for_time);
//|******************************************RESTART BEFORE CALENDERING*************************************************
    std::cout<<"Writing restart file ";
    simulator.write_restart_file(output_directory + "/pre_calendered_electrode_restart_file.res");
//|*********************************************************************************************************************
    //Calendaring process
    std::cout << "****************Initialize calendaring process ****************\n";
    bbox = simulator.get_bounding_box(); //get the XYZ max/min that contain all particles
    double h_3 = bbox[5]; //height of uppermost particle (Z-max)
    std::cout<<"Heigh of uppermoast particle: "<< h_3<< std::endl;
    if (h_3<1.2*calendering_height){
        h_3=1.2*calendering_height;
        std::cout<<"Heigh of uppermoast particle lower then 1.2 h_al: "<< h_3<< std::endl;
    }
    top_surface->move(-Vec3(0, 0,  calendering_height*3    - h_3-1.01*max_binder_thickness), Vec3(0, 0, 0)); //Move top surface to uppermost particle+binder thickness
    top_surface->set_velocity(Vec3(0,0,0.-surface_velocity/2));
    std::chrono::duration<double> compaction_time_2 {((h_3+1.01*max_binder_thickness-calendering_height) / (surface_velocity/2))};

    std::cout<<"Compaction time: "<< compaction_time_2.count()<< std::endl;
    run_for_time.reset(compaction_time_2);
    simulator.run(run_for_time);
    std::cout<<"Writing restart file ";
    simulator.write_restart_file(output_directory + "/compact_restart_file.res");
    std::cout << "****************Initialize unloading ****************\n";
    top_surface->set_velocity(Vec3(0, 0, surface_velocity/2));
//    mat->adhesive = true; //Enable adhesive when particles are compacted and before unloaded
    EngineType::SurfaceNormalForceLess zero_force(top_surface, 0.);
    // simulator.set_rotation(false);
    simulator.run(zero_force);
    top_surface->set_velocity(Vec3(0, 0, 0));
    auto top_surface_position_after_calendering = top_surface->get_points()[0].z();
    std::cout<<"Top surface position after calendering: "<< top_surface_position_after_calendering<< std::endl;
    std::cout << "****************Let rest for 1s ****************\n";
    run_for_time.reset(1e0s);
    simulator.run(run_for_time);
    bbox = simulator.get_bounding_box(); //get the XYZ max/min that contain all particles
    double Active_layer_height = bbox[5];
    std::cout<<"Height of active layer "<< Active_layer_height<< std::endl;
    //Calculate this from the position when the force on calendering surface is 0
    double NMC_porosity = 1-particle_volume/(box_side*box_side*top_surface_position_after_calendering);
//    double NMC_porosity = 1-particle_volume/(box_side*box_side*Active_layer_height) ;
    std::cout<<"NMC Porosity: "<< NMC_porosity<< std::endl;
    std::cout<<"Writing restart file ";
    simulator.write_restart_file(output_directory + "/calendered_electrode_restart_file.res");
    std::ofstream results_file;
    results_file.open(output_directory + "/calendered_electrode_results_file.dou");
    results_file << "Volume of particles=" << particle_volume << "\n";
    results_file<< "RVE side length= " << box_side << "\n";
    results_file << "NMC Porosity=" << NMC_porosity << "\n";
    results_file << "Height of active layer=" << Active_layer_height << "\n";
    results_file << "Top surface position=" << top_surface_position_after_calendering << "\n";
    results_file.close();
    */
}
