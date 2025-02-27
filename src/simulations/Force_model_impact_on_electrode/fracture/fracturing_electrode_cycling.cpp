//
// Adapted from electrode_cycling.cpp by Axel on 2025-02-19
//

#include "../../simulations.h"

#include <vector>

#include "../../../engine/engine.h"
#include "../../../contact_models/Positive_electrode/fracturing_swelling_elastic_plastic_binder_elastic_plastic_particle.h"
#include "../../../materials/electrode_material.h"


void DEM::fracturing_electrode_cycling(const std::string &settings_file_name)
{
    using namespace DEM;
    using ForceModel = fracturing_swelling_elastic_plastic_binder_elastic_plastic_particle;
    using ParticleType = FracturableSwellingSphericalParticle<ForceModel>;
    using EngineType = Engine<ForceModel, ParticleType>;
    using namespace std::chrono_literals;

    SimulationParameters parameters(settings_file_name);
    auto restart_file_name = parameters.get_parameter<std::string>("restart_file_name");
    auto output_directory = parameters.get_parameter<std::string>("output_dir");
    auto simulator = EngineType(restart_file_name);
    auto Calendering_output = simulator.get_output("output_0");
    simulator.remove_output(Calendering_output);

    auto mat = dynamic_cast<ElectrodeMaterial *>(simulator.get_material(0));
    mat->new_binder_contacts = false;

    std::chrono::duration<double> cycling_time{parameters.get_parameter<double>("cycling_time")};

    auto restart_array = parameters.get_vector<unsigned>("restart_array");

//==================SET TIME STEP AND MASS SCALING======================================================================
    double time_step = parameters.get_parameter<double>("time_step") * 1e-6;
    std::chrono::duration<double> time_step_us{time_step};
    std::cout << "Time step is:" << time_step_us.count() * 1E6 << " µs\n";
    simulator.set_time_increment(time_step_us);
    float mass_scaling_factor = parameters.get_parameter<float>("mass_scaling");
    std::cout << "Mass scaling is:" << mass_scaling_factor << "\n";
    simulator.set_mass_scale_factor(mass_scaling_factor);
//======================================================================================================================

//====================MAKE OUTPUT PRESCRIBED============================================================================
    double output_number = parameters.get_parameter<double>("output_number");
    std::chrono::duration<double> output_interval{cycling_time.count() / output_number};
    std::cout << "Output interval is: " << output_interval.count() << "s\n";
    std::cout << "Number of outputs are: " << cycling_time.count() / output_interval.count() << "s\n";
    auto cycling_output = simulator.create_output(output_directory, output_interval);
    cycling_output->print_particles = true;
    cycling_output->print_fractured_particles = true;
    cycling_output->print_surface_positions = true;
    cycling_output->print_kinetic_energy = true;
    cycling_output->print_contacts = true;
    cycling_output->print_surface_forces = true;
    cycling_output->print_fabric_force_tensor = true;
    cycling_output->print_periodic_bc = true;
    cycling_output->print_mirror_particles = true;
//======================================================================================================================

//======================Set swelling/scaling rate and time======================================================================

    auto swell_states = parameters.get_vector<double>("swell_states");
    auto material_scalings = parameters.get_vector<double>("material_scalings");

    if ( swell_states.size() != material_scalings.size() ) throw
        std::invalid_argument("swell_states and material_scaling arrays must be of same size!\n");
    std::cout << "Number of cycles for simulation: " <<  swell_states.size()/2 << "\n";

    // setup so swell_rate is constant (determined from first charge), charging time is determined from swell state and
    // material_scaling_rate is determined from scaling_factor and swellling_time.
    std::vector < std::tuple < int, double, std::chrono::duration < double>>> charging_direction_scale_rate_and_times;
    std::vector < std::tuple < double, double>> swell_and_scale_rates;
    double swell_state = 1.0;
    double material_scaling_state = 1.0;
//    int swell_direction = 0;
    for (int i = 0; i < swell_states.size(); ++i)
    {
        double to_be_swell_state = swell_states[i];
        double swell_state_diff = to_be_swell_state - swell_state;
        double swell_rate = swell_state_diff/ cycling_time.count();

        double to_be_material_scaling = material_scalings[i];
        double material_scaling_diff = to_be_material_scaling - material_scaling_state;
        double material_scaling_rate = material_scaling_diff/cycling_time.count();

        swell_and_scale_rates.push_back(std::make_tuple(swell_rate, material_scaling_rate));
        std::cout << "Charge step " << i << ", swell_rate = " << swell_rate <<
                  ", material_scaling_rate = " << material_scaling_rate << "\n";
        swell_state = to_be_swell_state;
        material_scaling_state = to_be_material_scaling;
    }

    //=MOVE TOP SURFACE FROM LAYER==================================================================================
    auto top_surface = simulator.get_surface<EngineType::PointSurfacePointer>("top_plate");
    //Remove calendering surface from RVE
    auto bbox = simulator.get_bounding_box(); //get the XYZ max/min that contain all particles
    double h_1 = bbox[5]; //height of uppermost particle (Z-max)
    double surface_removal_distance = h_1 + 1 - top_surface->get_points()[0].z();
    top_surface->move(Vec3(0, 0,  surface_removal_distance), Vec3(0, 0, 0));
    std::cout << "moving top surface with distance: " << surface_removal_distance << ", to height: " << h_1+1 << "\n";
    //==============================================================================================================

    //==REST LAYER BEFORE SIMULATION====================================================================================
    std::chrono::duration<double> resting_time{cycling_time.count() / 10.};
    std::cout << "Running for resting time: " << resting_time.count() << std::endl;
    EngineType::RunForTime run_for_resting_time(simulator, resting_time);
    simulator.run(run_for_resting_time);
    //==================================================================================================================

    //==ADD IN-PLANE STRAIN BOUNDARY CONDITION==========================================================================
    try
    {
        double in_plane_strain = parameters.get_parameter<double>("in_plane_strain");
        std::chrono::duration<double> in_plane_strain_loading_time{cycling_time.count() / 10.};
        double in_plane_strain_rate = in_plane_strain / in_plane_strain_loading_time.count();
        std::cout << "Applying an in-plane strain of:" << in_plane_strain << ", at a strain rate of:"
            << in_plane_strain_rate << ", for a time of:" << in_plane_strain_loading_time.count() <<"\n";
        simulator.set_periodic_boundary_condition_strain_rate('x',in_plane_strain_rate);
        auto deformable_surface = simulator.get_surface<EngineType::DeformablePointSurfacePointer>("bottom_plate");
        deformable_surface->set_in_plane_strain_rates(in_plane_strain_rate,0);

        EngineType::RunForTime run_for_periodic_BC_stretch(simulator, in_plane_strain_loading_time);
        simulator.run(run_for_periodic_BC_stretch);
        simulator.set_periodic_boundary_condition_strain_rate('x',0.);
        deformable_surface->set_in_plane_strain_rates(0.,0.);
    }
    catch (std::invalid_argument const &ex)
    {
        std::cout << ex.what()<<'\n';
        std::cout << "No in-plane strain found in input file.\n";
    }
    //==================================================================================================================

    //==ADD OUT OF PLANE STRESS BOUNDARY CONDITION======================================================================
    try
    {
        double out_of_plane_pressure = parameters.get_parameter<double>("out_of_plane_pressure");
        auto periodic_BCs = simulator.get_periodic_boundaries();
        auto RVE_out_of_plane_area = (periodic_BCs[0].max-periodic_BCs[0].min)*(periodic_BCs[1].max-periodic_BCs[1].min);
        double compaction_force = out_of_plane_pressure*RVE_out_of_plane_area;

        std::cout << "A out of plane pressure of: " << out_of_plane_pressure << ", and a out of plane area of: " <<
            RVE_out_of_plane_area << ", results in a compaction force of:" << compaction_force <<"\n";

//        double force_range = compaction_force/1e0;
//        double distance_range = 1e-1;
//        double time_range = 1e-3;//1e-2;
//        double velocity_range = distance_range / time_range;
//        double acceleration_range = velocity_range / time_range;

//        double surface_mass = parameters.get_parameter<double>("surface_mass");
//        double surface_mass = force_range / acceleration_range;
        double surface_mass {1e2};
//        double surface_damping_coefficient = parameters.get_parameter<double>("surface_damping_coefficient");
//        double surface_damping_coefficient = force_range / velocity_range;
        double surface_damping_coefficient {2.5E8}; // based on feeling

        std::cout << "The surface is given a mass of " << surface_mass <<
            ", and a damping coefficient of " << surface_damping_coefficient << "\n";

        // Move top surface to top of electrode layer
//        auto top_surface = simulator.get_surface<EngineType::PointSurfacePointer>("top_plate");
        auto bbox = simulator.get_bounding_box();
        double particle_height = bbox[5];

        double max_binder_thickness {0};
        for (auto& p: simulator.get_particles())
        {
            if (mat->binder_thickness_fraction*p->get_radius() > max_binder_thickness)
            {
                max_binder_thickness = mat->binder_thickness_fraction*p->get_radius();
            }
        }
        auto distance_to_move = top_surface->get_points()[0].z() - (particle_height + 1.01*max_binder_thickness) ;
        auto surface_velocity = particle_height / ( 5e0 * cycling_time.count() );
        top_surface->move(-Vec3(0,0, distance_to_move), -Vec3(0,0,surface_velocity));
        // Move top surface to that the desired compressive force is reached
        std::cout << "Compressing layer to a compaction force of: " << compaction_force
            << "with a velocity of : " << surface_velocity<< "\n";
        EngineType::SurfaceNormalForceGreater control_force(top_surface, compaction_force);
        simulator.run(control_force);
        top_surface->set_velocity(Vec3(0,0,0));

        if (!parameters.get_parameter<bool>("fixed_surface"))
        {

            std::cout << "Varying top layer with constant pressure.\n";
            // Add the prescribed force regulator
            top_surface->set_mass(surface_mass);
            top_surface->set_damping_coefficient(surface_damping_coefficient);
            auto amp_func = [compaction_force]() { return -compaction_force; };
            auto amp = std::make_shared<DEM::Amplitude>(amp_func);
            std::cout << "amp->value() = " << amp->value() << "\n";
            top_surface->set_force_amplitude(amp,'z');
//          simulator.set_force_control_on_surface(top_surface,'z'); //This function resets the surface force amplitude...
        }
        else std::cout << "Fixed top layer.\n";
    }
    catch (std::invalid_argument const &ex)
    {
        std::cout << ex.what() << '\n';
        std::cout << "No out of plane stress applied.\n";
    }

    //==================================================================================================================

    // Give array of after which cycles a restart should be written.
    unsigned counter{0};
    unsigned cycle{0};
    for (const auto &charge_iterator: swell_and_scale_rates)
    {
        std::cout << "Charge step " << counter << ", swell_rate = " << std::get<0>(charge_iterator) <<
                  ", material_scaling_rate = " << std::get<1>(charge_iterator) << "\n"<< std::endl;
        for (auto &p: simulator.get_particles())
        {
            p->set_swell_rate(std::get<0>(charge_iterator));
            p->set_material_scale_rate(std::get<1>(charge_iterator));
        }

        EngineType::RunForTime charging_run(simulator, cycling_time);
        simulator.run(charging_run);
        counter++;
        cycle = counter / 2;
        auto cycle_itr = std::find(restart_array.begin(), restart_array.end(), cycle);
        if ( cycle_itr != restart_array.end())
        {
            for (auto &p: simulator.get_particles())
            {
                p->set_swell_rate(0);
                p->set_material_scale_rate(0);
            }

            std::cout << "Writing restart file after cycle: " << cycle << "\n";
            restart_array.erase(cycle_itr);
            std::stringstream restart_file_name;
            restart_file_name << "/electrode_cycling_cycle_" << cycle << ".res";
            simulator.write_restart_file(output_directory + restart_file_name.str());
        }
    }
    std::cout << "**************** Simulation finalized ****************\n";
}