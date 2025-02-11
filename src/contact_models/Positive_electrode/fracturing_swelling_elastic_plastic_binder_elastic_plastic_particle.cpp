//
// Created by Axel on 2023-02-06.
// Adapted from elastic_plastic_binder_elastic_plastic_particle.cpp
//

#include "fracturing_swelling_elastic_plastic_binder_elastic_plastic_particle.h"

#include <cmath>
#include <iostream>
#include <random>
#include <vector>

#include "../../materials/electrode_material.h"

DEM::fracturing_swelling_elastic_plastic_binder_elastic_plastic_particle::fracturing_swelling_elastic_plastic_binder_elastic_plastic_particle(
        DEM::fracturing_swelling_elastic_plastic_binder_elastic_plastic_particle::ParticleType* particle1,
        DEM::fracturing_swelling_elastic_plastic_binder_elastic_plastic_particle::ParticleType* particle2,
        std::chrono::duration<double> dt)
{
    particle1_ = particle1;
    particle2_ = particle2;
    //extracting material data
    auto mat1 = dynamic_cast<const ElectrodeMaterial *>(particle1->get_material());
    auto mat2 = dynamic_cast<const ElectrodeMaterial *>(particle2->get_material());
    material = mat1;
    double R0 = 1. / (1. / particle1->get_radius() + 1. / particle2->get_radius());

    F_1_ = mat1->F_1_;
    alpha_1_ = mat1->alpha_1_;
    F_2_ = mat1->F_2_;
    alpha_2_ = mat1->alpha_2_;
    a_1_ = mat1->a_1_;
    beta_1_ = mat1->beta_1_;
    a_2_ = mat1->a_2_;
    beta_2_ = mat1->beta_2_;

    double E1 = mat1->E;
    v1 = mat1->nu;
    double vp1 = mat1->nup;
    double vp2 = mat2->nup;
    double Ep2 = mat2->Ep;
    double Ep1 = mat1->Ep;
    p1_fracture_strength_ = mat1->particle_fracture_strength_;
    p2_fracture_strength_ = mat2->particle_fracture_strength_;
    mu_particle_ = (mat1->mu + mat2->mu)/2;
    double Gp1 = Ep1/(2*(1+vp1));
    double Gp2 = Ep2/(2*(1+vp2));
//    double rhop = mat1->rhop;
    Ep_eff_ = 1./(((1-vp1*vp1)/Ep1)+((1-vp2*vp2)/Ep2));             //Effective particle youngs modulus
    br_ = mat1-> binder_radius_fraction * particle1->get_radius();
    bt_ = mat1-> binder_thickness_fraction * particle1->get_radius(); //Binder thickness determined by fraction of
    A = 3.1415 * br_ * br_;
    binder_yield_stress_ = mat1->binder_yield_stress_;
    particle_yield_stress_ = mat1->particle_yield_stress_;
    psi0_ = (1 - v1)/(1 + v1)/(1 - 2*v1)*E1*A/bt_; //The instantaneous value of the relaxation function for the binder
    psi0T_B_ = E1/bt_*A/2/(1+v1);
    kTp_ = 8/((2-vp1)/Gp1 + (2-vp2)/Gp2)*0.001*R0; //Tangential particle stiffness following Hertz contact for two particles, 0.001R0 represents
    adhesive_ = true;
//    new_binder_contacts_ = true;
    binder_contact_ = create_binder_contact(mat1);

    tau_i = mat1->tau_i;
    alpha_i = mat1->alpha_i;

    dt_ = dt.count(); //time increment
    M = mat1->M(); // size of Tau_i
    for (unsigned i=0; i!=M; ++i)
    {
    di_.push_back(0);
    ddi_.push_back(0);
    dti_Scalar.push_back(0);
    ddti_Scalar.push_back(0);
    ddti_.emplace_back(0., 0., 0.);
    dti_.emplace_back(0., 0., 0.);
    Ai.push_back(1-exp((-dt_/tau_i[i])));
    Bi.push_back((1-(tau_i[i]/dt_) *Ai[i]));
    }
}

DEM::fracturing_swelling_elastic_plastic_binder_elastic_plastic_particle::fracturing_swelling_elastic_plastic_binder_elastic_plastic_particle(
        DEM::fracturing_swelling_elastic_plastic_binder_elastic_plastic_particle::ParticleType* particle1,
        DEM::fracturing_swelling_elastic_plastic_binder_elastic_plastic_particle::SurfaceType* surface, std::chrono::duration<double> dt)
{
        particle1_ = particle1;
        auto mat1 = dynamic_cast<const ElectrodeMaterial *>(particle1->get_material());
        material = mat1;
        double R0 = particle1_->get_radius(); //Effective radius for one particle and one rigid surface

        F_1_ = mat1->F_1_;
        alpha_1_ = mat1->alpha_1_;
        F_2_ = mat1->F_2_;
        alpha_2_ = mat1->alpha_2_;
        a_1_ = mat1->a_1_;
        beta_1_ = mat1->beta_1_;
        a_2_ = mat1->a_2_;
        beta_2_ = mat1->beta_2_;

        double E1 = mat1->E;
        double v1 = mat1->nu;
        double vp1=mat1->nup;
        double Ep1 = mat1->Ep;
        p1_fracture_strength_ = mat1->particle_fracture_strength_;
        double Gp1 = Ep1/(2*(1+vp1));
        mu_particle_ = mat1->mu_wall;
//        double rhop = mat1->rhop;
        Ep_eff_ = 1./(((1-vp1*vp1)/Ep1)); //Effective Young's modulus for one particle and one rigid surface
        br_ = mat1-> binder_radius_fraction * particle1->get_radius();
        bt_ = mat1-> binder_thickness_fraction * particle1->get_radius();
        A = 3.1415 * br_ * br_;
        binder_yield_stress_ = mat1->binder_yield_stress_;
        particle_yield_stress_ = mat1->particle_yield_stress_;
        psi0_ = (1 - v1)/(1 + v1)/(1 - 2*v1)*E1*A/bt_; //The instantaneous value of the relaxation function for the binder
        psi0T_B_ = E1/bt_*A/2/(1+v1); //The instantaneous  value of the shear relaxation function for the binder
        kTp_ = 8/((2-vp1)/Gp1)*0.001*R0;                 //Tangential particle stiffness following Hertz contact for two particles, 0.001R0 represents
        adhesive_ = surface->adhesive();
//        new_binder_contacts_ = true;
        binder_contact_ = create_binder_contact(mat1);
        M = mat1->M(); // size of Tau_i
        tau_i = mat1->tau_i;
        alpha_i = mat1->alpha_i;
        dt_ = dt.count(); //time increment

        for (unsigned i=0; i!=M; ++i)
        {
            di_.push_back(0);
            ddi_.push_back(0);
            dti_Scalar.push_back(0);
            ddti_Scalar.push_back(0);
            ddti_.emplace_back(0., 0., 0.);
            dti_.emplace_back(0., 0., 0.);
            Ai.push_back(1-exp((-dt_/tau_i[i])));
            Bi.push_back((1-(tau_i[i]/dt_) *Ai[i]));
        }
}

DEM::fracturing_swelling_elastic_plastic_binder_elastic_plastic_particle::fracturing_swelling_elastic_plastic_binder_elastic_plastic_particle(
        DEM::fracturing_swelling_elastic_plastic_binder_elastic_plastic_particle::ParticleType* particle1,
        DEM::fracturing_swelling_elastic_plastic_binder_elastic_plastic_particle::ParticleType* particle2,
        std::chrono::duration<double>, const DEM::ParameterMap& parameters):
        particle1_(particle1),
        particle2_(particle2),
        psi0_(parameters.get_parameter<double>("psi0_")),
        psi0T_B_(parameters.get_parameter<double>("psi0T_B_")),
        F_1_(parameters.get_parameter<double>("F_1_")),
        alpha_1_(parameters.get_parameter<double>("alpha_1_")),
        F_2_(parameters.get_parameter<double>("F_2_")),
        alpha_2_(parameters.get_parameter<double>("alpha_2_")),
        a_1_(parameters.get_parameter<double>("a_1_")),
        beta_1_(parameters.get_parameter<double>("beta_1_")),
        a_2_(parameters.get_parameter<double>("a_2_")),
        beta_2_(parameters.get_parameter<double>("beta_2_")),
        p1_fracture_strength_(parameters.get_parameter<double>("p1_fracture_strength_")),
        p2_fracture_strength_(parameters.get_parameter<double>("p2_fracture_strength_")),
        fracture_degradation_factor_(parameters.get_parameter<double>("fracture_degradation_factor_")),
        kTp_(parameters.get_parameter<double>("kTp_")),
        bt_(parameters.get_parameter<double>("bt")),
        br_(parameters.get_parameter<double>("br")),
        A(parameters.get_parameter<double>("A")),
        Ep_eff_(parameters.get_parameter<double>("Ep_eff")),
        binder_yield_stress_(parameters.get_parameter<double>("binder_yield_stress_")),
        h_(parameters.get_parameter<double>("h")),
        hmax_(parameters.get_parameter<double>("hmax")),
        h_plast_(parameters.get_parameter<double>("h_plast_")),
        mu_particle_(parameters.get_parameter<double>("mu_particle")),
        particle_yield_stress_(parameters.get_parameter<double>("particle_yield_stress_")), //comment here
        bonded_(parameters.get_parameter<bool>("bonded_")),
        adhesive_(parameters.get_parameter<bool>("adhesive_")),
        binder_contact_(parameters.get_parameter<bool>("binder_contact")),
        binder_fracture_(parameters.get_parameter<bool>("binder_fracture")),
        particle_contact_(parameters.get_parameter<bool>("particle_contact")),
        a_0_(parameters.get_parameter<double>("a_0_")),
        dt_(parameters.get_parameter<double>("dt")),  // Time increment
        F_0_(parameters.get_parameter<double>("F_0_")),
        F_(parameters.get_parameter<double>("F")),
        dF_(parameters.get_parameter<double>("dF")),
        F_binder(parameters.get_parameter<double>("F_binder")),
        F_particle(parameters.get_parameter<double>("F_particle")),
        dFT_(parameters.get_vec3("dFT")),
        FT_(parameters.get_vec3("FT")),
        FT_binder_(parameters.get_vec3("FT_binder")),
        FT_binder_Scalar_(parameters.get_parameter<double>("FT_binder_Scalar_")),
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
        Bi.push_back(parameters.get_parameter<double>("B_" + std::to_string(i)));
        di_.push_back(parameters.get_parameter<double>("d_" + std::to_string(i)));
        ddi_.push_back(parameters.get_parameter<double>("dd_" + std::to_string(i)));
        dti_Scalar.push_back(parameters.get_parameter<double>("dti_Scalar"+ std::to_string(i)));
        ddti_Scalar.push_back(parameters.get_parameter<double>("ddti_Scalar"+ std::to_string(i)));
        dti_.push_back(parameters.get_vec3("dt_" + std::to_string(i)));
        ddti_.push_back(parameters.get_vec3("ddt_" + std::to_string(i)));
    }
}

DEM::fracturing_swelling_elastic_plastic_binder_elastic_plastic_particle::fracturing_swelling_elastic_plastic_binder_elastic_plastic_particle(
        DEM::fracturing_swelling_elastic_plastic_binder_elastic_plastic_particle::ParticleType* particle1,
        DEM::fracturing_swelling_elastic_plastic_binder_elastic_plastic_particle::SurfaceType * surface1,
        std::chrono::duration<double>, const DEM::ParameterMap& parameters):
        particle1_(particle1),
        psi0_(parameters.get_parameter<double>("psi0_")),
        psi0T_B_(parameters.get_parameter<double>("psi0T_B_")),
        F_1_(parameters.get_parameter<double>("F_1_")),
        alpha_1_(parameters.get_parameter<double>("alpha_1_")),
        F_2_(parameters.get_parameter<double>("F_2_")),
        alpha_2_(parameters.get_parameter<double>("alpha_2_")),
        a_1_(parameters.get_parameter<double>("a_1_")),
        beta_1_(parameters.get_parameter<double>("beta_1_")),
        a_2_(parameters.get_parameter<double>("a_2_")),
        beta_2_(parameters.get_parameter<double>("beta_2_")),
        p1_fracture_strength_(parameters.get_parameter<double>("p1_fracture_strength_")),
        fracture_degradation_factor_(parameters.get_parameter<double>("fracture_degradation_factor_")),
        kTp_(parameters.get_parameter<double>("kTp_")),
        Ep_eff_(parameters.get_parameter<double>("Ep_eff")),
        binder_yield_stress_(parameters.get_parameter<double>("binder_yield_stress_")),
        particle_yield_stress_(parameters.get_parameter<double>("particle_yield_stress_")),
        bt_(parameters.get_parameter<double>("bt")),
        br_(parameters.get_parameter<double>("br")),
        A(parameters.get_parameter<double>("A")),
        h_(parameters.get_parameter<double>("h")),
        hmax_(parameters.get_parameter<double>("hmax")),
        h_plast_(parameters.get_parameter<double>("h_plast_")),
        mu_particle_(parameters.get_parameter<double>("mu_particle")),
        bonded_(parameters.get_parameter<bool>("bonded_")),
        adhesive_(parameters.get_parameter<bool>("adhesive_")),
        binder_contact_(parameters.get_parameter<bool>("binder_contact")),
        binder_fracture_(parameters.get_parameter<bool>("binder_fracture")),
        particle_contact_(parameters.get_parameter<bool>("particle_contact")),
        a_0_(parameters.get_parameter<double>("a_0_")),
        dt_(parameters.get_parameter<double>("dt")),  // Time increment
        F_0_(parameters.get_parameter<double>("F_0_")),
        F_(parameters.get_parameter<double>("F")),
        dF_(parameters.get_parameter<double>("dF")),
        F_binder(parameters.get_parameter<double>("F_binder")),
        F_particle(parameters.get_parameter<double>("F_particle")),
        dFT_(parameters.get_vec3("dFT")),
        FT_(parameters.get_vec3("FT")),
        FT_binder_(parameters.get_vec3("FT_binder")),
        FT_binder_Scalar_(parameters.get_parameter<double>("FT_binder_Scalar_")),
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
        Bi.push_back(parameters.get_parameter<double>("B_" + std::to_string(i)));
        di_.push_back(parameters.get_parameter<double>("d_" + std::to_string(i)));
        ddi_.push_back(parameters.get_parameter<double>("dd_" + std::to_string(i)));
        dti_Scalar.push_back(parameters.get_parameter<double>("dti_Scalar"+ std::to_string(i)));
        ddti_Scalar.push_back(parameters.get_parameter<double>("ddti_Scalar"+ std::to_string(i)));
        dti_.push_back(parameters.get_vec3("dt_" + std::to_string(i)));
        ddti_.push_back(parameters.get_vec3("ddt_" + std::to_string(i)));
    }
}

void DEM::fracturing_swelling_elastic_plastic_binder_elastic_plastic_particle::update(
        double h, const DEM::Vec3& dt, const Vec3& drot, const DEM::Vec3& normal)
{
    rot_ += drot;
    F_ = update_normal_force(h);
    update_tangential_force(dt, normal);
    check_fracture();
}

unsigned DEM::fracturing_swelling_elastic_plastic_binder_elastic_plastic_particle::M;
const DEM::ElectrodeMaterial* DEM::fracturing_swelling_elastic_plastic_binder_elastic_plastic_particle::material;

double  DEM::fracturing_swelling_elastic_plastic_binder_elastic_plastic_particle::update_normal_force(double h)
{
    double dh = h-h_;
    h_ = h;
    if (h > hmax_)
    {
        hmax_ = h;
    }
// ==BINDER CONTACT MODEL===============================================================================================
    if(binder_contact_)
    {
        double h_def= hmax_ -  binder_yield_stress_*A/(psi0_*(1-alpha_i[0]));
                                                                      // Stiffness for determining the plastic overlap
                                                                      // is considering the relaxed binder stiffness
                                                                      // (i.e. considering the damping factor)
                                                                      // Plastic deformation of binder, where the binder
                                                                      // contact stiffness is constant in relation to
                                                                      // the initial binder thickness
        if (h_def < -bt_) h_def = -bt_;

        if (((h_ > h_def) && !particle_contact_) || (bonded_ && !particle_contact_))
        {
            double viscoelastic_summation = 0.;
            for (unsigned i = 0; i != M; i++)
            {
                ddi_[i] = Ai[i] * ((-h_def + h_) - di_[i]) + Bi[i] * dh;
                viscoelastic_summation += alpha_i[i] * ddi_[i];
                di_[i] += ddi_[i];
            }
            F_binder += psi0_ * (dh - viscoelastic_summation);

            double sigma_binder = F_binder/A;
            if (sigma_binder > binder_yield_stress_)
            {
                F_binder = binder_yield_stress_*A;
            }
            else if (sigma_binder < -binder_yield_stress_)
            {
                F_binder = -binder_yield_stress_*A;
            }

            //=EVALUATION OF FRACTURE IN TENSION WHICH REGARDS PLASTIC COMPRESSION AND LENGTH DECREASE OF BINDER========
            double binder_strain = (h_def - h_)/-h_def; // Evaluation of fracture in tension which regards
                                                        // the plastic deformation of the binder
//            double binder_strain = -(bt_ + h_)/bt_; // Evaluation of fracture in tension which does not regard
                                                    // the plastic deformation of the binder
            if ( binder_strain >= material->binder_fracture_strain_ and !binder_fracture_)
            {
                binder_fracture_ = true;
                bonded_ = false;
            }
            //==========================================================================================================
        }
        else if(particle_contact_) //Unloading of binder when particle contact is initiated
        {
            if (F_binder > 0 )
            {
                F_binder -= sqrt(pow((psi0_ * dh),2));
            }
            if (F_binder <= 0)
            {
                F_binder = 0;
            }
        }
        else //Removes any residual binder force if there is no binder contact
        {
            F_binder = 0;
            for(unsigned i = 0; i != M; i++)
            {
                di_[i] = 0;
            }
        }
        if ((F_binder > 0)  && adhesive() && !bonded_ && !binder_fracture_)
        {
            bonded_ = true;
        }
        else if (F_binder <= 0 && adhesive() && !bonded_ && !material->new_binder_contacts)
        {
            binder_contact_ = false;
        }
    }
    else
    {
        F_binder = 0;
        for(unsigned i = 0; i != M; i++)
        {
            di_[i] = 0;
        }
    }

// ==PARTICLE CONTACT MODEL=============================================================================================
    if (h_ > 0)
    {
//        if (material->bond_breaking)
//        {
//            particle_contact_ = true;
//            bonded_ = false;
//        }
        double R0 {0.};
        if (particle2_ != nullptr) R0 = 1./(1./particle1_->get_radius() + 1./particle2_->get_radius());
        else R0 = particle1_->get_radius();

        if ( particle1_->get_fracture()) fracture_degradation_factor_ = material->fracture_degradation_factor_;
        if ( particle2_ != nullptr)
        {
            if ( particle2_->get_fracture()) fracture_degradation_factor_ = material->fracture_degradation_factor_;
        }

        double mat_scaling {particle1_->get_material_scaling() * fracture_degradation_factor_};

        if (h_>=hmax_ || h_plast_ <= 0)
        {
//            hmax_ = h_;
            double h_norm = h_/R0;
            if (h_norm <= 0.25)
            {
                F_particle = pow(R0,2) * mat_scaling * particle_yield_stress_*
                        (F_1_ * pow(h_norm, alpha_1_) + F_2_ * pow(h_norm, alpha_2_));
                a_0_ = a_1_ * pow((h_norm), beta_1_) + a_2_ * pow((h_norm), beta_2_);
            }
            else if (h_norm > 0.25)
            {
                F_particle = pow(R0,2) * mat_scaling * particle_yield_stress_  *
                        ((F_1_ * pow((0.25), alpha_1_) + F_2_ * pow((0.25), alpha_2_)) +
                        (h_norm-0.25) * (F_1_ * alpha_1_ * pow((0.25), (alpha_1_ - 1)) +
                        F_2_ * alpha_2_ * pow((0.25), (alpha_2_ - 1))));
                a_0_ = (a_1_ * pow(0.25, beta_1_) + a_2_ * pow((0.25), beta_2_)) +
                        (h_norm-0.25) * (a_1_ * beta_1_ * pow((0.25), (beta_1_ - 1)) +
                        a_2_ * beta_2_ * pow((0.25), (beta_2_ - 1)));
            }
            h_plast_  = h_ - pow(pow(F_particle * (3./(4. * mat_scaling * Ep_eff_)),2)/R0,(1./3.));
        }
        // if overlap larger than elastic recover region initiation
        else if (h_ > h_plast_)
        {
            F_particle = 4./3. * mat_scaling * Ep_eff_ * pow(R0* pow((h_-h_plast_),3) ,1./2.);

            if(F_particle<0)
            {
                F_particle = 0;
            }
        }
        else {F_particle = 0;} // See how the hertz unloading should be performed
    }
    else
    {
        F_particle = 0;
    }

    if (adhesive() && bonded_ && !binder_fracture_)
    {
        return std::max(F_particle,0.)+F_binder;
    }
    else
    {
        return std::max(F_particle, 0.) + std::max(F_binder, 0.);
    }
}

void  DEM::fracturing_swelling_elastic_plastic_binder_elastic_plastic_particle::update_tangential_force(const DEM::Vec3& dt, const DEM::Vec3& normal)
{
    if (adhesive() && bonded_)
    {
        FT_binder_ -= dot_product(FT_binder_,normal)*normal;
        double uT_hist_mag = uT_.length();
        uT_ -= dot_product(uT_,normal)*normal;
        uT_ += dt;
        double uT_new_mag = uT_.length();
        double duT_mag = uT_new_mag-uT_hist_mag;
        double viscoelastic_sum = 0.0;
        for (unsigned i = 0; i !=M; i++)
        {
            ddti_Scalar[i] = Ai[i] * (uT_new_mag -dti_Scalar[i]) + Bi[i]*(duT_mag);
            viscoelastic_sum += alpha_i[i] * ddti_Scalar[i];
            dti_Scalar[i] += ddti_Scalar[i];
        }
        FT_binder_Scalar_ += psi0T_B_ * (duT_mag - viscoelastic_sum);
        if (uT_.length()==0 && FT_binder_.length()==0){
            FT_binder_.set_zero();
        }
        else if (uT_.length()==0)
        {
            FT_binder_ = FT_binder_Scalar_*-FT_.normal();
        }
        else
        {
            FT_binder_ = FT_binder_Scalar_*uT_.normal();
        }
    }

    else
    {
        FT_binder_Scalar_ = 0;
        uT_.set_zero();
        FT_binder_.set_zero();
        for (unsigned i = 0; i != M; ++i)
        {
            dti_[i].set_zero();
            dti_Scalar[i] = 0.0;
        }
    }
    FT_ = -FT_binder_;

    if (F_particle > 0.)
    {
        FT_part_ -= dot_product(FT_part_,normal)*normal;
        FT_part_ += kTp_*dt;
        if (FT_part_.length() > mu_particle_*F_particle)
        {
            FT_part_ = mu_particle_ * F_particle * FT_part_.normal();
        }
        else
        {
            FT_part_.set_zero();
        }
        FT_ -= FT_part_;
    }
}

void DEM::fracturing_swelling_elastic_plastic_binder_elastic_plastic_particle::check_fracture()
{
    if( !material->fracture) return;
    if( F_particle <= 0) return;
    if(!particle1_->get_fracture() and
    F_particle >= p1_fracture_strength_ * 2 * particle1_->get_radius() * 2 * particle1_->get_radius())
    {
        particle1_->set_fracture();
    }

    if( particle2_ != nullptr and
    !particle2_->get_fracture() and
    F_particle >= p2_fracture_strength_ * 2 * particle2_->get_radius() * 2 * particle2_->get_radius())
    {
        particle2_->set_fracture();
    }
}

std::string DEM::fracturing_swelling_elastic_plastic_binder_elastic_plastic_particle::restart_data() const {
        std::ostringstream ss;
        ss << named_print(dt_, "dt") << ", "
           << named_print(bt_, "bt") << ", "
           << named_print(br_, "br") << ", "
           << named_print(h_, "h") << ", "
           << named_print(hmax_, "hmax") << ", "
           << named_print(h_plast_, "h_plast_") << ", "
           << named_print(kTp_, "kTp_") << ", "
           << named_print(Ep_eff_, "Ep_eff") << ", "
           << named_print(binder_yield_stress_, "binder_yield_stress_") << ", "
           << named_print(particle_yield_stress_, "particle_yield_stress_") << ", "
           << named_print(F_, "F") << ", "
           << named_print(a_0_, "a_0_") << ", "
           << named_print(A, "A") << ", "
           << named_print(F_0_, "F_0_") << ", "
           << named_print(mu_particle_, "mu_particle") << ", "
           << named_print(dF_, "dF") << ", "
           << named_print(F_binder, "F_binder") << ", "
           << named_print(F_particle, "F_particle") << ", "
           << named_print(dFT_, "dFT") << ", "
           << named_print(FT_, "FT") << ", "
           << named_print(FT_binder_, "FT_binder") << ", "
           << named_print(FT_binder_Scalar_, "FT_binder_Scalar_") << ", "
           << named_print(FT_part_, "FT_part") << ", "
           << named_print(uT_, "uT") << ", "
           << named_print(rot_, "rot") << ", "
           << named_print(bonded_, "bonded_") << ", "
           << named_print(adhesive_, "adhesive_") << ", "
           << named_print(binder_contact_, "binder_contact") << ", "
           << named_print(binder_fracture_, "binder_fracture") << ", "
           << named_print(particle_contact_, "particle_contact") << ", "
           << named_print(fracture_degradation_factor_, "fracture_degradation_factor_") << ", "
           << named_print(p1_fracture_strength_, "p1_fracture_strength_") << ", "
           << named_print(p2_fracture_strength_, "p2_fracture_strength_") << ", "
           << named_print(psi0_, "psi0_") << ", "
           << named_print(psi0T_B_, "psi0T_B_") << ", "
           << named_print(F_1_, "F_1_") << ", "
           << named_print(alpha_1_, "alpha_1_") << ", "
           << named_print(F_2_, "F_2_") << ", "
           << named_print(alpha_2_, "alpha_2_") << ", "
           << named_print(a_1_, "a_1_") << ", "
           << named_print(beta_1_, "beta_1_") << ", "
           << named_print(a_2_, "a_2_") << ", "
           << named_print(beta_2_, "beta_2_") << ", "
           << named_print(M, "M");

        for (unsigned i=0; i != M; ++i) {
            ss <<  ", "
           << named_print(tau_i[i], "tau_" + std::to_string(i)) << ", "
           << named_print(alpha_i[i], "alpha_" + std::to_string(i)) << ", "
           << named_print(Ai[i], "A_" + std::to_string(i)) << ", "
           << named_print(Bi[i], "B_" + std::to_string(i)) << ", "
           << named_print(di_[i], "d_" + std::to_string(i)) << ", "
           << named_print(ddi_[i], "dd_" + std::to_string(i)) << ", "
           << named_print(dti_Scalar[i], "dti_Scalar" + std::to_string(i)) << ", "
           << named_print(ddti_Scalar[i], "ddti_Scalar" + std::to_string(i)) << ", "
           << named_print(dti_[i], "dt_" + std::to_string(i)) << ", "
           << named_print(ddti_[i], "ddt_" + std::to_string(i));
        }
        return ss.str();
    }

bool DEM::fracturing_swelling_elastic_plastic_binder_elastic_plastic_particle::create_binder_contact(const ElectrodeMaterial* mat)
{
    std::random_device random_device;
    std::default_random_engine rand_engine{ random_device() };
    std::uniform_real_distribution<double> distribution{0., 1.};
    double random_value = distribution(rand_engine);
    return random_value < mat->fraction_binder_contacts && mat->new_binder_contacts;
}

DEM::Vec3 DEM::fracturing_swelling_elastic_plastic_binder_elastic_plastic_particle::get_rolling_resistance_torque() const {
    //if (binder_contact_) {
    //    return -Rb_*Rb_*0.01*kB_*rot_;
    //}
    //else {
    return {0, 0, 0};
    //    }
}

std::string DEM::fracturing_swelling_elastic_plastic_binder_elastic_plastic_particle::get_output_string() const {
    std::stringstream ss;
    ss  << F_ << ", " << F_binder << ", " << F_particle << ", " << hmax_ << ", "
        << FT_.x() << ", "  << FT_.y() << ", " << FT_.z() << ", "
        << FT_binder_.x() << ", "  << FT_binder_.y() << ", " << FT_binder_.z() << ", "
        << FT_part_.x() << ", "  << FT_part_.y() << ", " << FT_part_.z() << ", "
        << binder_contact_ << ", " << bt_ << ", " << binder_fracture_;
    return ss.str();
}

void DEM::fracturing_swelling_elastic_plastic_binder_elastic_plastic_particle::set_increment(std::chrono::duration<double> dt) {
    dt_ = dt.count();
    Ai = {};
    Bi = {};
    for (unsigned i=0; i!=M; ++i) {
        Ai.push_back(1-exp((-dt_/tau_i[i])));
        Bi.push_back(tau_i[i]/dt_ *((dt_/tau_i[i])-Ai[i]));
    }
}

bool DEM::fracturing_swelling_elastic_plastic_binder_elastic_plastic_particle::adhesive() const
{
    return adhesive_ && material->adhesive;
}