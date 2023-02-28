//
// Created by Axel on 2021-12-06.
//
#include "../simulations.h"

#include <vector>

#include "../../engine/engine.h"
#include "../../contact_models/viscoelastic.h"
#include "../../contact_models/Binder_behavour_investigation/elastic_plastic_binder_rigid_plastic_particle.h"
#include "../../materials/electrode_material.h"
#include "../../utilities/file_reading_functions.h"
#include "../../utilities/filling_functions.h"

void DEM::electrode_natural_packing_rigid_perfect_plastic(const std::string& settings_file_name) {
    using namespace DEM;
    using ForceModel = elastic_plastic_binder_rigid_plastic_particle;
    using ParticleType = SphericalParticle<ForceModel>;
    using EngineType = Engine<ForceModel, ParticleType>;
    using namespace std::chrono_literals;

    SimulationParameters parameters(settings_file_name);

    float time_step = parameters.get_parameter<double>("time_step")*1e-6;
    std::chrono::duration<double> time_step_us {time_step};
    std::cout << "Time step is:" << time_step_us.count()*1E6 << "µs\n";
    EngineType simulator(time_step_us); //orig  1E0

    auto mass_scaling = parameters.get_parameter<double>("mass_scaling");
    simulator.set_mass_scale_factor(mass_scaling); //Orig 1E2
    std::cout << "Mass scaling is:" << mass_scaling << "\n";

    auto output_directory = parameters.get_parameter<std::string>("output_dir");
    auto N = parameters.get_parameter<double>("N"); //Number of particles
    auto particle_file = parameters.get_parameter<std::string>("radius_file");
    auto mat = simulator.create_material<ElectrodeMaterial>(4800);
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
    mat->particle_yield_stress_ = parameters.get_parameter<double>("particle_yield_stress_");
    mat->binder_radius_fraction = parameters.get_parameter<double>("binder_radius_fraction");
    mat->binder_thickness_fraction = parameters.get_parameter<double>("binder_thickness_fraction");
    mat->binder_yield_stress_ = parameters.get_parameter<double>("binder_yield_stress_");
    mat->binder_stiffness_coefficient = parameters.get_parameter<double>("binder_stiffness_coefficient");
    mat->fraction_binder_contacts = parameters.get_parameter<double>("fraction_binder_contacts");
    auto particle_density_at_filling = parameters.get_parameter<double>("filling_density");
    mat->active_particle_height = parameters.get_parameter<double>("active_particle_height");
//    mat->yield_displacement_coeff = parameters.get_parameter<double>("yield_displacement_coeff");

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
    std::cout << "box_height " << box_height << "\n";

    auto p1 = Vec3(-box_side / 2, -box_side / 2, 0);
    auto p2 = Vec3(box_side / 2, -box_side / 2, 0);
    auto p3 = Vec3(box_side / 2, box_side / 2, 0);
    auto p4 = Vec3(-box_side / 2, box_side / 2, 0);
    auto p5 = Vec3(-box_side / 2, -box_side / 2, box_height);
    auto p6 = Vec3(box_side / 2, -box_side / 2, box_height);
    auto p7 = Vec3(box_side / 2, box_side / 2, box_height);
    auto p8 = Vec3(-box_side / 2, box_side / 2, box_height);

    auto stiff_wall_fraction = 0.9;

    auto p1_stiff = Vec3(-stiff_wall_fraction*box_side / 2, -stiff_wall_fraction*box_side / 2, 0);
    auto p2_stiff = Vec3(stiff_wall_fraction*box_side / 2, -stiff_wall_fraction*box_side / 2, 0);
    auto p3_stiff = Vec3(stiff_wall_fraction*box_side / 2, stiff_wall_fraction*box_side / 2, 0);
    auto p4_stiff = Vec3(-stiff_wall_fraction*box_side / 2, stiff_wall_fraction*box_side / 2, 0);
    auto p5_stiff = Vec3(-stiff_wall_fraction*box_side / 2, -stiff_wall_fraction*box_side / 2, box_height);
    auto p6_stiff = Vec3(stiff_wall_fraction*box_side / 2, -stiff_wall_fraction*box_side / 2, box_height);
    auto p7_stiff = Vec3(stiff_wall_fraction*box_side / 2, stiff_wall_fraction*box_side / 2, box_height);
    auto p8_stiff = Vec3(-stiff_wall_fraction*box_side / 2, stiff_wall_fraction*box_side / 2, box_height);

    std::vector <Vec3> bottom_points{p1, p2, p3, p4};
    std::vector <Vec3> top_points{p8, p7, p6, p5};

    //side panels for non periodic BCs
    std::vector <Vec3> side_1{p2_stiff, p3_stiff, p7_stiff, p6_stiff};
    std::vector <Vec3> side_2{p3_stiff, p4_stiff, p8_stiff, p7_stiff};
    std::vector <Vec3> side_3{p4_stiff, p1_stiff, p5_stiff, p8_stiff};
    std::vector <Vec3> side_4{p1_stiff, p2_stiff, p6_stiff, p5_stiff};

    auto particle_positions = random_fill_box(-stiff_wall_fraction*box_side / 2, stiff_wall_fraction*box_side / 2, -stiff_wall_fraction*box_side / 2, stiff_wall_fraction*box_side / 2,
                                             0, 0+box_height, particle_radii, max_binder_thickness);

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

    for (std::size_t i = 0; i != particle_positions.size(); ++i) {
        simulator.create_particle(particle_radii[i], particle_positions[i], Vec3(0,0,0), mat);
    }

    mat->adhesive = false; //No adhesion of particles when initial packing
    double gravity = 1E1;
    simulator.set_gravity(Vec3(0, 0, -gravity)); //Use gravity for initial packing of particles
    std::cout << "max_binder_thickness: "<< max_binder_thickness <<"\n";
    simulator.setup(1.01*max_binder_thickness); //Size of box for detecting contacts between particles
//    simulator.set_rotation(false);

//============================CALCULATE FALL TIMES ==============================================================
    double fall_distance = box_height-mat->active_particle_height;
    std::cout << "Fall distance: "<< fall_distance<<" \n";
//    double  fall_time = pow(2*fall_distance/gravity ,0.5);
    std::chrono::duration<double> fall_time {pow(2*fall_distance/gravity ,0.5)};
    std::cout << "Fall time: "<< fall_time.count() <<" \n";
//============================Calculate output frequency================================================================
    int output_exp;
    if(std::log10(fall_time.count())>=0.0)
    {
        if(float(pow(10,(std::log10(fall_time.count())-int(std::log10(fall_time.count())))))>=5.0){output_exp = 1 +
                std::log10(fall_time.count());}else{output_exp = 0 + std::log10(fall_time.count());}
    }else{
        std::cout << "std::log10(fall_time.count()) " << std::log10(fall_time.count()) << "\n";
        if(float(pow(10,(std::log10(fall_time.count())-(int(std::log10(fall_time.count()))-1))))>=5.0){output_exp = 1 +
        std::log10(fall_time.count())-1;}else{output_exp = 0 + std::log10(fall_time.count())-1;}
    }
    std::cout << "fall time " << fall_time.count() << "\n";
    std::cout << "output exp " << output_exp << "\n";

   std::chrono::duration<double> output_interval {pow(10 ,output_exp-2)};
   std::cout << "Output interval: " << output_interval.count() << "\n";
//==================================OUTPUTS=============================================================================
    auto filling_output = simulator.create_output(output_directory , output_interval);
    filling_output->print_particles = true;
    filling_output->print_kinetic_energy = true;
    filling_output->print_surface_positions = true;
    filling_output->print_surface_forces = true;
    filling_output->print_contacts = true;
    filling_output->print_periodic_bc = true;
    filling_output->print_mirror_particles = true;
    filling_output->print_fabric_force_tensor=true;
//======================================================================================================================

    std::chrono::duration<double> acceleration_fall_time {fall_time.count()/10.0};
    std::cout << "****************Acceleration of particles**************** \n";
    std::cout << "Acceleration time: "<< acceleration_fall_time.count() <<" \n";
    EngineType::RunForTime Run_for_Particle_acceleration(simulator, acceleration_fall_time);
    simulator.run(Run_for_Particle_acceleration);

    simulator.set_gravity(Vec3(0,0,0));

    std::chrono::duration<double> fall_time_constant_vel {(fall_distance-gravity*pow(fall_time.count()/10.0,2)/1.0)/1.0/(gravity*fall_time.count()/10.0)};

    double pre_calendering_surface_velocity = 1 * gravity * fall_time.count()/10.0;
    std::cout << "Surface velocity: "<< pre_calendering_surface_velocity <<" \n";

    if (pre_calendering_surface_velocity*fall_time_constant_vel.count() >= box_height-3.2*mat->active_particle_height)
    {
    pre_calendering_surface_velocity = (box_height-3.2*mat->active_particle_height)/fall_time_constant_vel.count();
    std::cout << "Reducing surface velocity to: "<< pre_calendering_surface_velocity <<" \n";
    }

    top_surface->set_velocity(Vec3(0,0,-pre_calendering_surface_velocity));

    std::cout << "****************Falling without gravity**************** \n";
    std::cout << "No gravity fall time: "<< fall_time_constant_vel.count() <<" \n";
    EngineType::RunForTime Run_for_Particle_no_gravity(simulator, fall_time_constant_vel);
    simulator.run(Run_for_Particle_no_gravity);

    std::cout << "****************Initialize pre-calendering process**************** \n";
//  Move surface to uppermost particle
    simulator.set_gravity(Vec3(0,0,-gravity));
    auto bbox = simulator.get_bounding_box(); //get the XYZ max/min that contain all particles
    double h_1 = bbox[5]; //height of uppermost particle (Z-max)
    if (h_1<3.1*mat->active_particle_height){
        h_1=3.1*mat->active_particle_height;
        std::cout<<"Height of uppermost particle lower then 3.1 h_al: "<< h_1<< std::endl;
    }
    else{
        std::cout<<"Height of uppermost particle: "<< h_1<< std::endl;
    }

    top_surface->move(-Vec3(0, 0,  top_surface->get_points()[0][2] - h_1-1.01*max_binder_thickness), Vec3(0, 0, 0)); //Move top surface to uppermost partile+binder thickness

    top_surface->set_velocity(Vec3(0,0,-2 * pre_calendering_surface_velocity));

    std::chrono::duration<double> compaction_time_pre_cal {((h_1+1.01*max_binder_thickness - mat->active_particle_height*3) / (2 * pre_calendering_surface_velocity))};

    std::cout<<"Pre-calendering time: "<< compaction_time_pre_cal.count()<< std::endl;
    EngineType::RunForTime Run_for_Pre_calendering_time(simulator, compaction_time_pre_cal);

    simulator.run(Run_for_Pre_calendering_time);

    top_surface->set_velocity(Vec3(0,0,0));

    EngineType::RunForTime Run_for_rest_time(simulator, fall_time);
    std::cout<<"Resting for: "<< fall_time.count()<< std::endl;
    simulator.run(Run_for_rest_time);

    std::cout << "****************Wall removal**************** \n";
    // Stop all the particles
    for (auto& p: simulator.get_particles())
    {
        p->set_velocity(Vec3(0,0,0));
    }
//=======================================PERIODIC BC:S===============================================================
    simulator.add_periodic_boundary_condition('x', -box_side/2, box_side/2);
    simulator.add_periodic_boundary_condition('y', -box_side/2, box_side/2);
//=====================================================================================================================

// =====================MOVE THE STIFF SURFACE TO INITIATE THE PERIODIC BC:S ==========================================
    side1_surface->move(Vec3(5*box_side,0,0), Vec3(0, 0, 0));
    side2_surface->move(Vec3(0,5*box_side,0), Vec3(0,0,0));
    side3_surface->move(-Vec3(5*box_side,0,0), Vec3(0,0,0));
    side4_surface->move(-Vec3(0,5*box_side,0), Vec3(0,0,0));
//=====================================================================================================================

    EngineType::RunForTime Run_for_initiation_of_periodic_BCs(simulator,fall_time);
    std::cout << "Running for: "<< (fall_time).count() <<" \n";
    simulator.run(Run_for_initiation_of_periodic_BCs);

    // Stop all the particles
    for (auto& p: simulator.get_particles())
    {
        p->set_velocity(Vec3(0,0,0));
    }

    Run_for_initiation_of_periodic_BCs.reset(fall_time);
    std::cout << "Stopping particles and running for: "<< (fall_time).count() <<" \n";
    simulator.run(Run_for_initiation_of_periodic_BCs);

    std::cout << "****************Adhesive on and stopping particles**************** \n";
    // Stop all the particles
    for (auto& p: simulator.get_particles())
    {
        p->set_velocity(Vec3(0,0,0));
    }
    mat->adhesive = true; // Activate adhesion before calendering starts

    EngineType::RunForTime Run_for_adhesive_resting_time(simulator, fall_time);
    std::cout << "Running for: "<< (fall_time).count() <<" \n";
    simulator.run(Run_for_adhesive_resting_time);

    //turn off gravity
    simulator.set_gravity(Vec3(0, 0, 0));
    std::cout << "****************Turning of gravity**************** \n";

    EngineType::RunForTime Run_for_gravity_removal_resting_time(simulator,fall_time);
    std::cout << "Running for: "<< (fall_time).count() <<" \n";
    simulator.run(Run_for_gravity_removal_resting_time);
//========================================== WRITE RESTART BEFORE CALENDERING===========================================
    std::cout<<"Writing restart file ";
    simulator.write_restart_file(output_directory + "/pre_calendered_electrode_restart_file.res");
//======================================================================================================================
}
