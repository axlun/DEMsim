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

        double E{ 0. };    // E0 for the binder material
        double nu{ 0. };   // nu for the binder material
        double binder_yield_stress_{ 1.0e99 };  //Yield strength for the binder material
        double particle_yield_stress_{ 1.0e99}; //Yield strength for the particle
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
    };
}
#endif