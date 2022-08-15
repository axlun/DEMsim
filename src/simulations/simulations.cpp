//
// Created by erolsson on 2018-09-18.
//

#include <map>

#include "simulations.h"

std::map<std::string, DEM::SimulationFunctionPtr> DEM::valid_simulations() {
    // Add all valid simulation routines here
    return std::map<std::string, SimulationFunctionPtr>{
            {"gyratory_compaction",                     DEM::gyratory_compaction},
            {"closed_die_compaction",                   DEM::closed_die_compaction},
            {"cube_die_compaction",                     DEM::cube_die_compaction},
            {"cathode_compaction",                      DEM::cathode_compaction},
            {"contact_tester",                          DEM::contact_tester},
            {"contact_tester_elastic_perfect_plastic_material", DEM::contact_tester_elastic_perfect_plastic_material},
            {"binder_contact_tester_viscoelastic_binder_El_Pl_particles", DEM::binder_contact_tester_viscoelastic_binder_El_Pl_particles},
            {"particle_contact_tester_viscoelastic_binder_El_Pl_particles", DEM::particle_contact_tester_viscoelastic_binder_El_Pl_particles},
            {"wall_contact_tester_viscoelastic_binder_El_Pl_particles", DEM::wall_contact_tester_viscoelastic_binder_El_Pl_particles},
            {"contact_tester_viscoelastic_binder_El_Pl_particles", DEM::contact_tester_viscoelastic_binder_El_Pl_particles},
            {"cyclic_triaxial",                         DEM::cyclic_triaxial},
            {"proctor_test",                            DEM::proctor_test},
            {"stone_compaction",                        DEM::stone_compaction},
            {"electrode_compaction",                    DEM::electrode_compaction},
            {"electrode_mechanical_test",               DEM::electrode_mechanical_test},
            {"electrode_cylinder_compaction",           DEM::electrode_cylinder_compaction},
            {"periodic_bc_tester",                      DEM::periodic_bc_tester},
            {"periodic_bc_simulation",                  DEM::periodic_bc_simulation},
            {"filling_periodic_box",                    DEM::filling_periodic_box},
            {"Cathode_mechanical_simulations",          DEM::Cathode_mechanical_simulations},
            {"battery_rve_compaction",                  DEM::battery_rve_compaction},
            {"restart_electrode",                       DEM::restart_electrode},
            {"porous_electrode_rve",                    DEM::porous_electrode_rve},
            {"electrode_calendering",                   DEM::electrode_calendering},
            {"electrode_calendering_hertz",             DEM::electrode_calendering_hertz},
            {"restart_test",                            DEM::restart_test},
            {"electrode_mechanical_loading",            DEM::electrode_mechanical_loading},
            {"electrode_mechanical_loading_hertz",      DEM::electrode_mechanical_loading_hertz},
            {"restart_electrode_calendering_hertz",     DEM::restart_electrode_calendering_hertz},
            {"restart_electrode_calendering",           DEM::restart_electrode_calendering}
    };
}

