//
// Created by elahe on 2019-11-14.
//
#include "material_base.h"

#include <vector>

#ifndef DEMSIM_VISCOELASTICMATERIAL_H
#define DEMSIM_VISCOELASTICMATERIAL_H

namespace DEM {
    class ParameterMap;
    class ElectrodeMaterial : public MaterialBase {
    public:
        ElectrodeMaterial(unsigned id_number, double density) : MaterialBase(id_number, density) {}
        explicit ElectrodeMaterial(const ParameterMap& parameters);
        ~ElectrodeMaterial() override = default;
        [[nodiscard]] std::string restart_data() const override;
        double F_1_{0.};
        double alpha_1_{0.};
        double F_2_{0.};
        double alpha_2_{0.};
        double a_1_{0.};
        double beta_1_{0.};
        double a_2_{0.};
        double beta_2_{0.};
        double E{ 0. };    // E0 for the binder material
        double nu{ 0. };   // nu for the binder material
        double binder_yield_stress_{ 1.0e99 };  //Yield strength for the binder material
        double binder_fracture_strain_ { 1.0e99}; //Fracture strain for binder
        double particle_yield_stress_{ 1.0e99}; //Yield strength for the particle
        double particle_fracture_strength_{ 1.0e99};
        double fracture_degradation_factor_{1.0};
        double Ep{ 0. };   // Young's modulus for the particles
        double nup{ 0. };  // Poisson's ratio for the binder material
        double rhop{ 0. }; //density of particles

        double bt { 0. };   // Binder thickness

        // double active_particle_height

        double yield_displacement_coeff { 8.59e-3 }; // Yield displacement of particle

        double N=0.;


        // double contact;
        std::vector<double> tau_i {};
        std::vector<double> alpha_i {};
        double fraction_binder_contacts{ 0. };
        double binder_radius_fraction{ 0. };
        double binder_thickness_fraction{ 0. };
        double binder_stiffness_coefficient{ 1. };

        double kT;
        double mu{ 0. };
        double mu_wall{ 0. };
        double mu_binder { 0. };

        [[nodiscard]] unsigned M() const { return tau_i.size(); }

        double active_particle_height=0.;
//        bool bond_breaking = true;
        bool adhesive = true;
        bool fracture = true;
        bool new_binder_contacts = true;
    };
}
#endif