//
// Created by erolsson on 12/08/2020.
//

#include "electrode_material.h"
#include <iostream>

#include "../utilities/file_reading_functions.h"
#include "../utilities/printing_functions.h"

DEM::ElectrodeMaterial::ElectrodeMaterial(const ParameterMap& parameters) :
    DEM::MaterialBase(parameters),
    E(parameters.get_parameter<double>("E")),
    nu(parameters.get_parameter<double>("nu")),
    binder_yield_stress_(parameters.get_parameter<double>("binder_yield_stress_")),//Binder yield strength
    Ep(parameters.get_parameter<double>("Ep")),
    nup(parameters.get_parameter<double>("nup")),
    rhop(parameters.get_parameter<double>("rhop")),
    bt(parameters.get_parameter<double>("bt")),
    particle_yield_stress_(parameters.get_parameter<double>("particle_yield_stress_")),
    yield_displacement_coeff(parameters.get_parameter<double>("yield_coeff")),
    tau_i(),
    alpha_i(),
    fraction_binder_contacts(parameters.get_parameter<double>("fraction_binder_contacts")),
    binder_radius_fraction(parameters.get_parameter<double>("binder_radius_fraction")),
    binder_thickness_fraction(parameters.get_parameter<double>("binder_thickness_fraction")),
    binder_stiffness_coefficient(parameters.get_parameter<double>("binder_stiffness_coefficient")),

    kT(parameters.get_parameter<double>("kT")),
    mu(parameters.get_parameter<double>("mu")),
    mu_wall(parameters.get_parameter<double>("mu_wall")),
//    bond_breaking(parameters.get_parameter<double>("bond_breaking")),
    adhesive(parameters.get_parameter<bool>("adhesive"))
{
    auto M = parameters.get_parameter<std::size_t>("M");
    for (std::size_t i = 0; i != M; ++i) {
        tau_i.push_back(parameters.get_parameter<double>("tau_" + std::to_string(i)));
        alpha_i.push_back(parameters.get_parameter<double>("alpha_" + std::to_string(i)));
    }
}

std::string DEM::ElectrodeMaterial::restart_data() const {
    std::ostringstream ss;
    ss << named_print("electrode_material", "type") << ", "
       << MaterialBase::restart_data() << ", "
       << named_print(E, "E") << ", "
       << named_print(nu, "nu") << ", "
       << named_print(binder_yield_stress_, "binder_yield_stress_") << ", "
       << named_print(nup, "nup") << ", "
       << named_print(Ep, "Ep") << ", "
       << named_print(rhop, "rhop") << ", "
       << named_print(fraction_binder_contacts,"fraction_binder_contacts")<< ", "
       << named_print(binder_thickness_fraction, "binder_thickness_fraction")<< ", "
       << named_print(binder_stiffness_coefficient, "binder_stiffness_coefficient")<< ", "
       << named_print(yield_displacement_coeff, "yield_coeff") << ", "
       << named_print(binder_radius_fraction, "binder_radius_fraction") << ", "
       << named_print(particle_yield_stress_, "particle_yield_stress_") << ", "
       << named_print(bt, "bt") << ", "
       << named_print(kT, "kT") << ", "
       << named_print(mu, "mu") << ", "
       << named_print(mu_wall, "mu_wall") << ", "
       << named_print(M(), "M")<< ", "
        << named_print(bond_breaking, "bond_breaking")<< ", "
       << named_print(adhesive, "adhesive") << ", ";

    for (std::size_t i = 0; i != M(); ++i) {
        ss << ", " << named_print(tau_i[i], "tau_"+ std::to_string(i)) << ", "
           << named_print(alpha_i[i], "alpha_"+ std::to_string(i));
    }
    return ss.str();
}

