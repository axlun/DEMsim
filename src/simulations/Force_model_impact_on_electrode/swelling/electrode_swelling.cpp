//
// Adapted from electrode_resting_el_pl_binder_el_pl_particle.cpp on 2024-06-26.
//

#include "../../simulations.h"

#include <vector>

#include "../../../engine/engine.h"
#include "../../../contact_models/Positive_electrode/swelling_elastic_plastic_binder_elastic_plastic_particle.h"
#include "../../../materials/electrode_material.h"


void DEM::electrode_swelling(const std::string &settings_file_name)
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
//    mat->new_binder_contacts = false;

    std::chrono::duration<double> swelling_time{parameters.get_parameter<double>("swelling_time")};

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
    double output_number =parameters.get_parameter<double>("output_number");
    std::chrono::duration<double> output_interval {swelling_time.count()/output_number};
    std::cout << "Output interval is: " << output_interval.count() << "s\n";
    std::cout << "Number of outputs are: " << swelling_time.count() / output_interval.count() << "s\n";
    auto swelling_output = simulator.create_output(output_directory, output_interval);
    swelling_output->print_particles = true;
    swelling_output->print_surface_positions = true;
    swelling_output->print_kinetic_energy = true;
    swelling_output->print_contacts = true;
    swelling_output->print_surface_forces = true;
    swelling_output->print_fabric_force_tensor = true;
    swelling_output->print_periodic_bc = true;
    swelling_output->print_mirror_particles = true;
//======================================================================================================================

//======================Set swelling/scaling rate and time======================================================================

    double swell_state {1};
    auto swell_states = parameters.get_vector<double>("swell_states");
    double swell_rate = (swell_states.back() - swell_state) / swelling_time.count();
    std::cout << "swell_rate = " << swell_rate << std::endl;

    auto material_scaling = parameters.get_parameter<double>("material_scaling");
    double material_scaling_rate = (material_scaling - 1) / swelling_time.count();

    std::vector <std::chrono::duration<double>> swell_times;
    for ( const auto &ss: swell_states)
    {
        if (swell_rate == 0.0)
        {
            swell_times.push_back(swelling_time);
            std::cout << "No swell rate!" << std::endl;
            break;
        }
        std::chrono::duration<double> st {(ss - swell_state) / swell_rate};
        swell_times.push_back(st);
        swell_state = ss;
    }
    std::chrono::duration<double> resting_time{swelling_time.count()/10.};
    std::cout << "Running for resting time: " << resting_time.count() << std::endl;
    EngineType::RunForTime run_for_resting_time(simulator, resting_time);
    simulator.run(run_for_resting_time);
    // for (auto &p: simulator.get_particles())
    // {
    //     p->set_swell_rate(swell_rate);
    // }
    //  Set swell rate for particles
    //  loop for all swelling times  and run for swelling time, stop swell rate, write restart and ad swell rate
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