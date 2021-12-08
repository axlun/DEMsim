//
// Created by Axel on 2021-11-17.
//

#include "viscoelastic_binder_El_Pl_particles.h"

#include <cmath>
#include <iostream>
#include <random>
#include <vector>

#include "../../materials/electrode_material.h"

DEM::viscoelastic_binder_El_Pl_particles::viscoelastic_binder_El_Pl_particles(DEM::viscoelastic_binder_El_Pl_particles::ParticleType* particle1,
                                         DEM::viscoelastic_binder_El_Pl_particles::ParticleType* particle2,
                                         std::chrono::duration<double> dt)
 {
    //extracting material data
    auto mat1 = dynamic_cast<const ElectrodeMaterial *>(particle1->get_material());
    auto mat2 = dynamic_cast<const ElectrodeMaterial *>(particle2->get_material());
    material = mat1;

    R0_ = 1. / (1. / particle1->get_radius() + 1. / particle2->get_radius());

    double E1 = mat1->E;
    double v1 = mat1->nu;
    double E2 = mat1->E;
    double v2 = mat1->nu;
    double vp1 = mat1->nup;
    double vp2 = mat2->nup;
    double Ep2 = mat2->Ep;
    double Ep1 = mat1->Ep;
    double Ep_eff = 1./(((1-vp1*vp1)/Ep1)+((1-vp2*vp2)/Ep2));
    double br_ = mat1-> binder_radius_fraction * particle1->get_radius();
    bt_ = mat1-> binder_thickness_fraction * particle1->get_radius(); //Binder thickness determined by fraction of
//  std::cout << "Binder thickness:" << bt_ << std::endl;                      //particle radius
    double A = DEM::pi*br_*br_;
    kb_coeff = mat1->binder_stiffness_coefficient;
    psi0_ = kb_coeff*(1 - v1)/(1 + v1)/(1 - 2*v1)*E1*A/bt_; //The instantaneous value of the relaxation function for the binder
    kp_ = (4./3)*Ep_eff* sqrt(R0_); //Particle stiffness following Hertz contact for two particles
    yield_h_ = mat1->yield_displacement_coeff*R0_;


    binder_contact_ = create_binder_contact(mat1);
    adhesive_ = true;



    tau_i = mat1->tau_i;
    alpha_i = mat1->alpha_i;

    dt_ = dt.count(); //time increment
    M = mat1->M(); // size of Tau_i
//  std::cout << "M:" << M << std::endl;
    for (unsigned i=0; i!=M; ++i)
    {
    di_.push_back(0);
    ddi_.push_back(0);
    Ai.push_back(1-exp((-dt_/tau_i[i])));
    Bi.push_back((1-(tau_i[i]/dt_) *Ai[i]));
    }
}

DEM::viscoelastic_binder_El_Pl_particles::viscoelastic_binder_El_Pl_particles(
        DEM::viscoelastic_binder_El_Pl_particles::ParticleType* particle1,
        DEM::viscoelastic_binder_El_Pl_particles::SurfaceType* surface, std::chrono::duration<double> dt)
{
        auto mat1 = dynamic_cast<const ElectrodeMaterial *>(particle1->get_material());
        material = mat1;
        R0_ = particle1->get_radius(); //Effective radius for one particle and one rigid surface
        double E1 = mat1->E;
        double v1 = mat1->nu;
        double vp1=mat1->nup;
        double Ep1 = mat1->Ep;
        double Ep_eff = 1./(((1-vp1*vp1)/Ep1)); //Effective Young's modulus for one particle and one rigid surface
        double br_ = mat1-> binder_radius_fraction * particle1->get_radius();
        bt_ = mat1-> binder_thickness_fraction * particle1->get_radius();
        double A = DEM::pi*br_*br_;
        psi0_ = (1 - v1)/(1 + v1)/(1 - 2*v1)*E1*A/bt_; //The instantaneous value of the relaxation function for the binder
        kp_ = (4./3)*Ep_eff* sqrt(R0_); //Particle stiffness following Hertz contact for two particles
        yield_h_ = mat1->yield_displacement_coeff*R0_;
        adhesive_ = surface->adhesive();
        binder_contact_ = create_binder_contact(mat1);


        M = mat1->M(); // size of Tau_i

        tau_i = mat1->tau_i;
        alpha_i = mat1->alpha_i;

        dt_ = dt.count(); //time increment

        for (unsigned i=0; i!=M; ++i)
        {
            di_.push_back(0);
            ddi_.push_back(0);
            Ai.push_back(1-exp((-dt_/tau_i[i])));
            Bi.push_back((1-(tau_i[i]/dt_) *Ai[i]));
        }
}


DEM::viscoelastic_binder_El_Pl_particles::viscoelastic_binder_El_Pl_particles(
        DEM::viscoelastic_binder_El_Pl_particles::ParticleType* particle1,
        DEM::viscoelastic_binder_El_Pl_particles::ParticleType* particle2,
        std::chrono::duration<double>, const DEM::ParameterMap& parameters):
        psi0_(parameters.get_parameter<double>("psi0_")),
        kp_(parameters.get_parameter<double>("kparticle")),
        R0_(parameters.get_parameter<double>("R0")),
        //Rb_(parameters.get_parameter<double>("Rb")),
        bt_(parameters.get_parameter<double>("bt")),
        h_(parameters.get_parameter<double>("h")),
        yield_h_(parameters.get_parameter<double>("yield_h")),
        hmax_(parameters.get_parameter<double>("hmax")),
        mu_particle_(parameters.get_parameter<double>("mu_particle")),
        adhesive_(parameters.get_parameter<bool>("adhesive_")),
        binder_contact_(parameters.get_parameter<bool>("binder_contact")),
        fractured_(parameters.get_parameter<bool>("fractured")),
        dt_(parameters.get_parameter<double>("dt")),  // Time increment
        F_(parameters.get_parameter<double>("F")),
        dF_(parameters.get_parameter<double>("dF")),
        F_binder(parameters.get_parameter<double>("F_binder")),
        F_particle(parameters.get_parameter<double>("F_particle")),
        dFT_(parameters.get_vec3("dFT")),
        FT_(parameters.get_vec3("FT")),
        FT_binder_(parameters.get_vec3("FT_binder")),
        FT_part_(parameters.get_vec3("FT_part")),
        uT_(parameters.get_vec3("uT")),
        rot_(parameters.get_vec3("rot"))
        {
            material = dynamic_cast<const ElectrodeMaterial *>(particle1->get_material());
            M = parameters.get_parameter<unsigned>("M");
            for (unsigned i = 0; i != M; ++i) {
                tau_i.push_back(parameters.get_parameter<double>("tau_" + std::to_string(i)));
                alpha_i.push_back(parameters.get_parameter<double>("alpha_" + std::to_string(i)));
                Ai.push_back(parameters.get_parameter<double>("A_" + std::to_string(i)));
                Bi.push_back(parameters.get_parameter<double>("A_" + std::to_string(i)));
                di_.push_back(parameters.get_parameter<double>("d_" + std::to_string(i)));
                ddi_.push_back(parameters.get_parameter<double>("dd_" + std::to_string(i)));
                dti_.push_back(parameters.get_vec3("dt_" + std::to_string(i)));
                ddti_.push_back(parameters.get_vec3("ddt_" + std::to_string(i)));
            }
        }

DEM::viscoelastic_binder_El_Pl_particles::viscoelastic_binder_El_Pl_particles(
        DEM::viscoelastic_binder_El_Pl_particles::ParticleType* particle1,
        DEM::viscoelastic_binder_El_Pl_particles::SurfaceType * surface1,
        std::chrono::duration<double>, const DEM::ParameterMap& parameters):
        psi0_(parameters.get_parameter<double>("psi0_")),
        kp_(parameters.get_parameter<double>("kparticle")),
        R0_(parameters.get_parameter<double>("R0")),
        //Rb_(parameters.get_parameter<double>("Rb")),
        bt_(parameters.get_parameter<double>("bt")),
        h_(parameters.get_parameter<double>("h")),
        yield_h_(parameters.get_parameter<double>("yield_h")),
        hmax_(parameters.get_parameter<double>("hmax")),
        mu_particle_(parameters.get_parameter<double>("mu_particle")),
        adhesive_(parameters.get_parameter<bool>("adhesive_")),
        binder_contact_(parameters.get_parameter<bool>("binder_contact")),
        fractured_(parameters.get_parameter<bool>("fractured")),
        dt_(parameters.get_parameter<double>("dt")),  // Time increment
        F_(parameters.get_parameter<double>("F")),
        dF_(parameters.get_parameter<double>("dF")),
        F_binder(parameters.get_parameter<double>("F_binder")),
        F_particle(parameters.get_parameter<double>("F_particle")),
        dFT_(parameters.get_vec3("dFT")),
        FT_(parameters.get_vec3("FT")),
        FT_binder_(parameters.get_vec3("FT_binder")),
        FT_part_(parameters.get_vec3("FT_part")),
        uT_(parameters.get_vec3("uT")),
        rot_(parameters.get_vec3("rot"))
{
    material = dynamic_cast<const ElectrodeMaterial *>(particle1->get_material());
    M = parameters.get_parameter<unsigned>("M");
    for (unsigned i = 0; i != M; ++i) {
        tau_i.push_back(parameters.get_parameter<double>("tau_" + std::to_string(i)));
        alpha_i.push_back(parameters.get_parameter<double>("alpha_" + std::to_string(i)));
        Ai.push_back(parameters.get_parameter<double>("A_" + std::to_string(i)));
        Bi.push_back(parameters.get_parameter<double>("A_" + std::to_string(i)));
        di_.push_back(parameters.get_parameter<double>("d_" + std::to_string(i)));
        ddi_.push_back(parameters.get_parameter<double>("dd_" + std::to_string(i)));
        dti_.push_back(parameters.get_vec3("dt_" + std::to_string(i)));
        ddti_.push_back(parameters.get_vec3("ddt_" + std::to_string(i)));
    }
}

void DEM::viscoelastic_binder_El_Pl_particles::update(double h, const DEM::Vec3& dt, const Vec3& drot, const DEM::Vec3& normal)
{
//    std::cout << "New iteration" << std::endl;
    F_ = update_normal_force(h);
//    std::cout << "F_:"<< F_ << std::endl;
//    std::cout << "F_binder:"<< F_binder << std::endl;
//    std::cout << "F_particle:"<< F_particle << std::endl;
//    std::cout << "binder_contact_:"<< binder_contact_ << std::endl;
//    std::cout << "bonded_:"<< bonded_ << std::endl;
//    std::cout << "adhesive_:"<< adhesive_ << std::endl;
}

unsigned DEM::viscoelastic_binder_El_Pl_particles::M;
const DEM::ElectrodeMaterial* DEM::viscoelastic_binder_El_Pl_particles::material;

double  DEM::viscoelastic_binder_El_Pl_particles::update_normal_force(double h)
{
   // std::cout << "*****New Iteration*****" << h_ << std::endl;

    double dh = h-h_;
    h_ = h;
//    std::cout << "h:" << h_ << std::endl;
//    std::cout << "dh:" << dh << std::endl; //change

    if (h > hmax_) {
        hmax_ = h;
    }

    if(binder_contact_)
    {
        if((h>-bt_)||bonded_)
        {

        double viscoelastic_summation = 0.;
        for(unsigned i = 0; i != M; i++)
        {
            ddi_[i] = Ai[i] * ((bt_+h_) - di_[i]) + Bi[i] * dh;
//            std::cout << "Ai:" << Ai[i] << std::endl;
//            std::cout << "Bi:" << Bi[i] << std::endl;
//            std::cout << "di_:" << di_[i] << std::endl;
//            std::cout << "ddi_:" << ddi_[i] << std::endl;
            viscoelastic_summation += alpha_i[i]*ddi_[i];
            di_[i] += ddi_[i];
        }
//        std::cout << "Viscoelastic_summation:" << -viscoelastic_summation << std::endl;
//        std::cout << "Elastic force increment:" << psi0_*dh << std::endl;

        F_binder += psi0_*(dh-viscoelastic_summation);
//        std::cout << "F_binder:" << F_binder << std::endl;
        }
        else
        {
            F_binder = 0;
            for(unsigned i = 0; i != M; i++)
            {
                di_[i] = 0;
            }
        }
        if (F_binder > 0  && adhesive())
        {
            bonded_ = true;
        }
    }
    if (h_ > 0) {
        if (h > yield_h_ && h >= hmax_)
        {
            F_particle += 1.5 * kp_ * sqrt(yield_h_) * dh; //Tangent of Hertz at yield displacement?
        }
        else {
            F_particle += 1.5*kp_*sqrt(h_)*dh;
        }
    }
    else
    {
        F_particle = 0;
    }
    if (adhesive() && bonded_)
    {
        return std::max(F_particle,0.)+F_binder;
    }
    else {
        return std::max(F_particle, 0.) + std::max(F_binder, 0.);
    }
}

std::string DEM::viscoelastic_binder_El_Pl_particles::restart_data() const {
        std::ostringstream ss;
        ss << named_print(dt_, "dt") << ", "
           << named_print(bt_, "bt") << ", "
           << named_print(h_, "h") << ", "
           << named_print(hmax_, "hmax") << ", "
           << named_print(yield_h_, "yield_h") << ", "
           //<< named_print(k_, "k") << ", "
           << named_print(kp_, "kp_") << ", "
           << named_print(R0_, "R0") << ", "
           //<< named_print(Rb_, "Rb") << ", "
           << named_print(F_, "F") << ", "
           << named_print(mu_particle_, "mu_particle") << ", "
           //<< named_print(mu_binder_, "mu_binder") << ", "
           << named_print(dF_, "dF") << ", "
           << named_print(F_binder, "F_binder") << ", "
           << named_print(F_particle, "F_particle") << ", "
           << named_print(dFT_, "dFT") << ", "
           << named_print(FT_, "FT") << ", "
           << named_print(FT_binder_, "FT_binder") << ", "
           << named_print(FT_part_, "FT_part") << ", "
           << named_print(uT_, "uT") << ", "
           << named_print(rot_, "rot") << ", "
           << named_print(bonded_, "activated") << ", "
           << named_print(adhesive_, "adhesive_") << ", "
           << named_print(binder_contact_, "binder_contact") << ", "
           << named_print(fractured_, "fractured") << ", "
           << named_print(psi0_, "psi0_") << ", "
           << named_print(kT_B_, "kT_B_") << ", "
           << named_print(M, "M");

        for (unsigned i=0; i != M; ++i) {
            ss <<  ", "
               << named_print(tau_i[i], "tau_" + std::to_string(i)) << ", "
               << named_print(alpha_i[i], "alpha_" + std::to_string(i)) << ", "
               << named_print(Ai[i], "A_" + std::to_string(i)) << ", "
               << named_print(Bi[i], "B_" + std::to_string(i)) << ", "
               << named_print(di_[i], "d_" + std::to_string(i)) << ", "
               << named_print(ddi_[i], "dd_" + std::to_string(i)) << ", "
               << named_print(dti_[i], "dt_" + std::to_string(i)) << ", "
               << named_print(ddti_[i], "ddt_" + std::to_string(i));
        }
        return ss.str();
    }





bool DEM::viscoelastic_binder_El_Pl_particles::create_binder_contact(const ElectrodeMaterial* mat)
{
    std::random_device random_device;
    std::default_random_engine rand_engine { random_device() };
    std::uniform_real_distribution<double> distribution{0., 1.};
    double random_value = distribution(rand_engine);
    if (random_value < mat->fraction_binder_contacts)
    {
        return true;
    }
    return false;
}

DEM::Vec3 DEM::viscoelastic_binder_El_Pl_particles::get_rolling_resistance_torque() const {
    //if (binder_contact_) {
    //    return -Rb_*Rb_*0.01*kB_*rot_;
    //}
    //else {
    return {0, 0, 0};
    //    }
}

std::string DEM::viscoelastic_binder_El_Pl_particles::get_output_string() const {
    std::stringstream ss;
    ss  << F_ << ", " << F_binder << ", " << F_particle << ", " << hmax_ << ", "
        << FT_.x() << ", "  << FT_.y() << ", " << FT_.z() << ", "
        << FT_binder_.x() << ", "  << FT_binder_.y() << ", " << FT_binder_.z() << ", "
        << FT_part_.x() << ", "  << FT_part_.y() << ", " << FT_part_.z() << ", "
        << binder_contact_ << ", " << bt_;
    return ss.str();
}

void DEM::viscoelastic_binder_El_Pl_particles::set_increment(std::chrono::duration<double> dt) {
    dt_ = dt.count();
    Ai = {};
    Bi = {};
    for (unsigned i=0; i!=M; ++i) {
        Ai.push_back(1-exp((-dt_/tau_i[i])));
        Bi.push_back(tau_i[i]/dt_ *((dt_/tau_i[i])-Ai[i]));
    }
}

bool DEM::viscoelastic_binder_El_Pl_particles::adhesive() const
{
    return adhesive_ && material->adhesive;
}