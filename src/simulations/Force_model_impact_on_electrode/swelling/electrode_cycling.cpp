//
// Created by Axel on 2024-08-19
//

#include "../../simulations.h"

#include <vector>

#include "../../../engine/engine.h"
#include "../../../contact_models/Positive_electrode/swelling_elastic_plastic_binder_elastic_plastic_particle.h"
#include "../../../materials/electrode_material.h"


void DEM::electrode_cycling(const std::string &settings_file_name)
{
    using namespace DEM;
    using ForceModel = swelling_elastic_plastic_binder_elastic_plastic_particle;
    using ParticleType = SwellingSphericalParticle<ForceModel>;
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

//    auto compression_pressure = parameters.get_parameter<double>("compression_pressure");
//    auto surface_adhesion = parameters.get_parameter<bool>("surface_adhesion");
//    auto surface_velocity = parameters.get_parameter<double>("surface_velocity");

//==================SET TIME STEP AND MASS SCALING=======================================================================
    double time_step = parameters.get_parameter<double>("time_step") * 1e-6;
    std::chrono::duration<double> time_step_us{time_step};
    std::cout << "Time step is:" << time_step_us.count() * 1E6 << " µs\n";
    simulator.set_time_increment(time_step_us);
    float mass_scaling_factor = parameters.get_parameter<float>("mass_scaling");
    std::cout << "Mass scaling is:" << mass_scaling_factor << "\n";
    simulator.set_mass_scale_factor(mass_scaling_factor);
//======================================================================================================================

//====================MAKE OUTPUT PRESCRIBED==============================================================================
    double output_number = parameters.get_parameter<double>("output_number");
    std::chrono::duration<double> output_interval{cycling_time.count() / output_number};
    std::cout << "Output interval is: " << output_interval.count() << "s\n";
    std::cout << "Number of outputs are: " << cycling_time.count() / output_interval.count() << "s\n";
    auto cycling_output = simulator.create_output(output_directory, output_interval);
    cycling_output->print_particles = true;
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

    if ( swell_states.size() != material_scalings.size() ) throw std::invalid_argument("swell_states and material_scaling arrays must be of same size!\n");
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

    std::chrono::duration<double> resting_time{cycling_time.count() / 10.};
    std::cout << "Running for resting time: " << resting_time.count() << std::endl;
    EngineType::RunForTime run_for_resting_time(simulator, resting_time);
//    simulator.run(run_for_resting_time);

    //==PLACE TOP SURFACE ON ACTIVE LAYER===============================================================================
    //move top_surface "closer" to active layer
    /*
    auto top_surface = simulator.get_surface<EngineType::PointSurfacePointer>("top_plate");

    auto surface_area = pow(2*top_surface->get_points()[0][0], 2);
    std::cout << "surface_area = " << surface_area << std::endl;

    auto compression_force = compression_pressure * 1E6 * surface_area;
    std::cout << "compression_force = " << compression_force << std::endl;

    auto bottom_surface = simulator.get_surface<EngineType::DeformablePointSurfacePointer>("bottom_plate");

    top_surface->set_velocity(Vec3(0,0,-surface_velocity));
    top_surface->set_adhesive(surface_adhesion);
    std::cout << "**************** Moving top surface with velocity "<< surface_velocity << ", to pressure " <<
        compression_pressure <<" ****************\n";

    EngineType::SurfaceNormalForceGreater compaction_force(top_surface, compression_force);
    simulator.run(compaction_force);

    top_surface->set_velocity(Vec3(0, 0, 0));
    */
    //==================================================================================================================
    //=MOVE TOP SURFACE FROM LAYER======================================================================================
    auto calendering_surface = simulator.get_surface<EngineType::PointSurfacePointer>("top_plate");
    //Remove calendering surface from RVE
    auto bbox = simulator.get_bounding_box(); //get the XYZ max/min that contain all particles
    double h_1 = bbox[5]; //height of uppermost particle (Z-max)
    double surface_removal_distance = h_1 + 1 - calendering_surface->get_points()[0].z();
    calendering_surface->move(Vec3(0, 0,  surface_removal_distance), Vec3(0, 0, 0));
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