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

    std::chrono::duration<double> cycling_time{parameters.get_parameter<double>("cycling_time")};

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

    auto swell_state = parameters.get_parameter<double>("swell_state");
    double swell_rate = (swell_state - 1) / cycling_time.count();
    std::cout << "swell_rate = " << swell_rate << std::endl;

    auto material_scaling = parameters.get_parameter<double>("material_scaling");
    double material_scaling_rate = (material_scaling - 1) / cycling_time.count();

    auto states_of_charge = parameters.get_vector<double>("states_of_charge");

    // Set up so swelling/scaling rate is constant, and the time betweeel 0 - 100% SOC is prescribed.
    std::vector < std::pair < int, std::chrono::duration < double>>> charging_directions_and_times;
    double SOC = 0.0;
    int charge_direction = 0;
    for (int i = 0; i < states_of_charge.size(); ++i)
    {
        double to_be_SOC = states_of_charge[i];
        double SOC_diff = to_be_SOC - SOC;
        charge_direction = -1;
        if (SOC_diff >= 0) charge_direction = 1;
        std::chrono::duration<double> charge_time{cycling_time * abs(SOC_diff) / 100.0};
        charging_directions_and_times.push_back(std::make_pair(charge_direction, charge_time));
        SOC = to_be_SOC;
    }

    std::chrono::duration<double> resting_time{cycling_time.count() / 10.};
    std::cout << "Running for resting time: " << resting_time.count() << std::endl;
    EngineType::RunForTime run_for_resting_time(simulator, resting_time);
    simulator.run(run_for_resting_time);

    //==PLACE TOP SURFACE ON ACTIVE LAYER===============================================================================
    //move top_surface "closer" to active layer
    auto top_surface = simulator.get_surface<EngineType::PointSurfacePointer>("top_plate");

    auto surface_area = pow(2*top_surface->get_points()[0][0], 2);
    std::cout << "surface_area = " << surface_area << std::endl;
    auto compression_pressure = parameters.get_parameter<double>("compression_pressure");

    auto compression_force = compression_pressure * 1E6 * surface_area;
    std::cout << "compression_force = " << compression_force << std::endl;

    auto bottom_surface = simulator.get_surface<EngineType::DeformablePointSurfacePointer>("bottom_plate");

    auto surface_velocity = parameters.get_parameter<double>("surface_velocity");
    top_surface->set_velocity(Vec3(0,0,-surface_velocity));
    auto surface_adhesion = parameters.get_parameter<bool>("surface_adhesion");
    top_surface->set_adhesive(surface_adhesion);
    std::cout << "**************** Moving top surface with velocity "<< surface_velocity << ", to pressure " <<
        compression_pressure <<" ****************\n";

    EngineType::SurfaceNormalForceGreater compaction_force(top_surface, compression_force);
    simulator.run(compaction_force);

    top_surface->set_velocity(Vec3(0, 0, 0));
    //==================================================================================================================

    unsigned counter{0};
    for (const auto &charge_iterator: charging_directions_and_times)
    {
        std::cout << "Charge step " << counter << ", charging direction = " << charge_iterator.first <<
                  ", charging time = " << charge_iterator.second.count() << std::endl;
        for (auto &p: simulator.get_particles())
        {
            p->set_swell_rate(charge_iterator.first * swell_rate);
            p->set_material_scale_rate(charge_iterator.first * material_scaling_rate);
        }

        EngineType::RunForTime charging_run(simulator, charge_iterator.second);
        simulator.run(charging_run);
        counter++;

    }
    std::cout << "**************** Simulation finalized ****************\n";
}
/*

    for (int i = 0; i < swell_times.size(); i++)
    {
        std::cout << "Swelling step " << i << ", swelling for time: " << swell_times[i].count() << std::endl;
        for (auto &p: simulator.get_particles())
        {
            p->set_swell_rate(swell_rate);
            p->set_material_scale_rate(material_scaling_rate);
        }
        EngineType::RunForTime run_for_swelling_time(simulator, swell_times[i]);
        simulator.run(run_for_swelling_time);
        std::cout << "Swelling finished, resting for " << resting_time.count() << "s" << std::endl;
        for (auto &p: simulator.get_particles())
        {
            p->set_swell_rate(0.);
            p->set_material_scale_rate(0.);
        }
        run_for_resting_time.reset(resting_time);
        simulator.run(run_for_resting_time);
        std::cout << "**************** Simulation finished ****************" << std::endl;
        std::cout << "Writing restart file " << std::endl;
        std::stringstream restart_file_name;
        restart_file_name << "/electrode_swell_state_" << swell_states[i] << ".res";
        simulator.write_restart_file(output_directory + restart_file_name.str());
    }
    std::cout << "**************** Simulation finalized ****************\n";
}
 */