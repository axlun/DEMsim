//
// Created by erolsson on 2018-09-18.
//

#ifndef SIMULATIONS_H
#define SIMULATIONS_H

#include <map>
#include <string>

namespace DEM {
    using SimulationFunctionPtr = void (*)(const std::string&);

    void gyratory_compaction(const std::string& settings_file_name);
    void closed_die_compaction(const std::string& settings_file_name);
    void cathode_compaction(const std::string& settings_file_name);
    void cube_die_compaction(const std::string& settings_file_name);
    void material_scaling_particle_tester(const std::string& settings_file_name);
    void fracturing_particle_tester(const std::string& settings_file_name);
    void fracturing_particle_plate_compression(const std::string& settings_file_name);
    void fracturing_particle_cube_die_compaction(const std::string& settings_file_name);
    void fracturing_particle_periodic_compaction(const std::string& settings_file_name);
    void fracturing_periodic_packing(const std::string& settings_file_name);
    void fracturing_electrode_calendering(const std::string& settings_file_name);
    void fracturing_electrode_mechanical_loading(const std::string& settings_file_name);
    void swelling_particle_tester(const std::string& settings_file_name);
    void swelling_particle_wall_tester(const std::string& settings_file_name);
    void swelling_cube_die_compaction(const std::string& settings_file_name);
    void material_scaling_cube_die_compaction(const std::string& settings_file_name);
    void material_scaling_periodic_compaction(const std::string& settings_file_name);
    void cube_die_packing_and_swelling_wc(const std::string& settings_file_name);
    void swelling_cube_die(const std::string& settings_file_name);
    void swelling_periodic_compaction(const std::string& settings_file_name);
    void swelling_periodic_packing(const std::string& settings_file_name);
    void swelling_electrode_calendering(const std::string& settings_file_name);
    void swelling_electrode_mechanical_loading(const std::string& settings_file_name);
    void electrode_cycling(const std::string& settings_file_name);
    void periodic_bc_restart_tester(const std::string& settings_file_name);
    void contact_tester(const std::string& settings_file_name);
    void contact_tester_elastic_perfect_plastic_material(const std::string& settings_file_name);
    void contact_tester_viscoelastic_binder_El_Pl_particles(const std::string& settings_file_name);
    void binder_tangential_contact_tester_elastic_plastic_binder_hertz_particles(const std::string& settings_file_name);
    void particle_normal_contact_tester_elastic_plastic_binder_rigid_plastic_particles(const std::string& settings_file_name);
    void particle_normal_contact_tester_elastic_plastic_binder_elastic_plastic_particles(const std::string& settings_file_name);
    void particle_force_overlap_el_pl_particles(const std::string& settings_file_name);
    void particle_force_overlap_hertz_particles(const std::string& settings_file_name);
    void particle_to_wall_drop_contact_tester_elastic_plastic_binder_elastic_plastic_particles(const std::string& settings_file_name);
    void binder_wall_tangential_contact_tester_elastic_plastic_binder_hertz_particles(const std::string& settings_file_name);
    void normal_contact_tester_elastic_plastic_binder_hertz_particles(const std::string& settings_file_name);
    void particle_to_wall_contact_tester_elastic_plastic_binder_hertz_particles(const std::string& settings_file_name);
    void particle_to_wall_drop_contact_tester_elastic_plastic_binder_rigid_plastic_particles(const std::string& settings_file_name);
    void binder_contact_tester_viscoelastic_binder_El_Pl_particles(const std::string& settings_file_name);
    void particle_contact_tester_viscoelastic_binder_El_Pl_particles(const std::string& settings_file_name);
    void wall_contact_tester_viscoelastic_binder_El_Pl_particles(const std::string& settings_file_name);
    void cyclic_triaxial(const std::string& settings_file_name);
    void proctor_test(const std::string& settings_file_name);
    void stone_compaction(const std::string& settings_file_name);
    void electrode_compaction(const std::string& settings_file_name);
    void electrode_swelling(const std::string& settings_file_name);
    void electrode_cylinder_compaction(const std::string& settings_file_name);
    void periodic_bc_tester(const std::string& settings_file_name);
    void periodic_bc_simulation(const std::string& settings_file_name);
    void periodic_bc_restart_tester(const std::string& settings_file_name);
    void filling_periodic_box(const std::string& settings_file_name);
    void battery_rve_compaction(const std::string &settings_file_name);
    void Cathode_mechanical_simulations(const std::string& settings_file_name);
    void restart_electrode(const std::string& settings_file_name);
    void porous_electrode_rve(const std::string& settings_file_name);
    void natural_packing_hertz(const std::string& settings_file_name);
    void electrode_calendering_hertz(const std::string& settings_file_name);
    void electrode_resting_hertz(const std::string& settings_file_name);
    void electrode_natural_packing_hertz(const std::string& settings_file_name);
    void electrode_natural_packing_el_pl_binder_el_pl_particle(const std::string& settings_file_name);
    void electrode_natural_periodic_packing_el_pl_binder_el_pl_particle(const std::string& settings_file_name);
    void electrode_natural_packing_rigid_perfect_plastic(const std::string& settings_file_name);
    void electrode_calendering(const std::string& settings_file_name);
    void restart_test(const std::string& settings_file_name);
    void restart_electrode_calendering(const std::string& settings_file_name);
    void restart_file_tester(const std::string& settings_file_name);
    void restart_electrode_calendering_hertz(const std::string& settings_file_name);
    void restart_electrode_calendering_el_pl_binder_el_pl_particle(const std::string& settings_file_name);
    void restart_electrode_calendering_rigid_plastic(const std::string& settings_file_name);
    void electrode_mechanical_loading(const std::string& settings_file_name);
    void electrode_mechanical_loading_hertz(const std::string& settings_file_name);
    void electrode_mechanical_loading_rigid_plastic(const std::string& settings_file_name);
    void electrode_mechanical_loading_el_pl_binder_el_pl_particle(const std::string& settings_file_name);
    void electrode_relaxation_el_pl_binder_el_pl_particle(const std::string& settings_file_name);
    void electrode_resting_el_pl_binder_el_pl_particle(const std::string& settings_file_name);
    std::map<std::string, SimulationFunctionPtr> valid_simulations();


    void electrode_mechanical_test(const std::string &settings_file_name);

}

#endif //SIMULATIONS_H
