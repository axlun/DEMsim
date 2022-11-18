//
// Created by Axel on 2021-11-17.
//

#include "elastic_plastic_binder_hertz_plastic_particle.h"

#include <cmath>
#include <iostream>
#include <random>
#include <vector>

#include "../../materials/electrode_material.h"

DEM::elastic_plastic_binder_hertz_plastic_particle::elastic_plastic_binder_hertz_plastic_particle(DEM::elastic_plastic_binder_hertz_plastic_particle::ParticleType* particle1,
                                                                                                                  DEM::elastic_plastic_binder_hertz_plastic_particle::ParticleType* particle2,
                                                                                                                  std::chrono::duration<double> dt)
 {
    //extracting material data
    auto mat1 = dynamic_cast<const ElectrodeMaterial *>(particle1->get_material());
    auto mat2 = dynamic_cast<const ElectrodeMaterial *>(particle2->get_material());
    material = mat1;
    R0_ = 1. / (1. / particle1->get_radius() + 1. / particle2->get_radius());
    double E1 = mat1->E;
    v1 = mat1->nu;
    double vp1 = mat1->nup;
    double vp2 = mat2->nup;
    double Ep2 = mat2->Ep;
    double Ep1 = mat1->Ep;
    particle_yield_stress_ = mat1->particle_yield_stress_;
    mu_particle_ = (mat1->mu + mat2->mu)/2;
    double Gp1 = Ep1/(2*(1+vp1));
    double Gp2 = Ep2/(2*(1+vp2));
//    double rhop = mat1->rhop;
    Ep_eff_ = 1./(((1-vp1*vp1)/Ep1)+((1-vp2*vp2)/Ep2));             //Effective particle youngs modulus
    br_ = mat1-> binder_radius_fraction * particle1->get_radius();
    bt_ = mat1-> binder_thickness_fraction * particle1->get_radius(); //Binder thickness determined by fraction of
    A = DEM::pi*br_*br_;
    kb_coeff = mat1->binder_stiffness_coefficient; //reduction of binder stiffness due to geometry
    binder_yield_stress_ = mat1->binder_yield_stress_;
    yield_h_ = mat1->yield_displacement_coeff*R0_; //Yield of particles

    psi0_ = kb_coeff*(1 - v1)/(1 + v1)/(1 - 2*v1)*E1*A/bt_; //The instantaneous value of the relaxation function for the binder
    psi0T_B_ = E1/bt_*A/2/(1+v1);
    kTp_ = 8/((2-vp1)/Gp1 + (2-vp2)/Gp2)*0.001*R0_; //Tangential particle stiffness following Hertz contact for two particles, 0.001R0 represents
    kp_ = (4./3)*Ep_eff_* sqrt(R0_); //Particle stiffness following Hertz contact for two particles
    adhesive_ = true;
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

DEM::elastic_plastic_binder_hertz_plastic_particle::elastic_plastic_binder_hertz_plastic_particle(
        DEM::elastic_plastic_binder_hertz_plastic_particle::ParticleType* particle1,
        DEM::elastic_plastic_binder_hertz_plastic_particle::SurfaceType* surface, std::chrono::duration<double> dt)
{
        auto mat1 = dynamic_cast<const ElectrodeMaterial *>(particle1->get_material());
        material = mat1;
        R0_ = particle1->get_radius(); //Effective radius for one particle and one rigid surface
        double E1 = mat1->E;
        double v1 = mat1->nu;
        double vp1=mat1->nup;
        double Ep1 = mat1->Ep;
        double Gp1 = Ep1/(2*(1+vp1));
        particle_yield_stress_ = mat1->particle_yield_stress_;
        mu_particle_ = mat1->mu_wall;
//        double rhop = mat1->rhop;
        Ep_eff_ = 1./(((1-vp1*vp1)/Ep1)); //Effective Young's modulus for one particle and one rigid surface
        br_ = mat1-> binder_radius_fraction * particle1->get_radius();
        bt_ = mat1-> binder_thickness_fraction * particle1->get_radius();
        A = DEM::pi*br_*br_;
        kb_coeff = mat1->binder_stiffness_coefficient; //reduction of binder stiffness due to geometry
        binder_yield_stress_ = mat1->binder_yield_stress_;
        yield_h_ = mat1->yield_displacement_coeff*R0_; //Yield of particles
        psi0_ = kb_coeff*(1 - v1)/(1 + v1)/(1 - 2*v1)*E1*A/bt_; //The instantaneous value of the relaxation function for the binder
        psi0T_B_ = E1/bt_*A/2/(1+v1); //The instantaneous  value of the shear relaxation function for the binder
        kTp_ = 8/((2-vp1)/Gp1)*0.001*R0_;                 //Tangential particle stiffness following Hertz contact for two particles, 0.001R0 represents
        kp_ = (4./3)*Ep_eff_* sqrt(R0_); //Particle stiffness following Hertz contact for two particles
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
            dti_Scalar.push_back(0);
            ddti_Scalar.push_back(0);
            ddti_.emplace_back(0., 0., 0.);
            dti_.emplace_back(0., 0., 0.);
            Ai.push_back(1-exp((-dt_/tau_i[i])));
            Bi.push_back((1-(tau_i[i]/dt_) *Ai[i]));
        }
}

DEM::elastic_plastic_binder_hertz_plastic_particle::elastic_plastic_binder_hertz_plastic_particle(
        DEM::elastic_plastic_binder_hertz_plastic_particle::ParticleType* particle1,
        DEM::elastic_plastic_binder_hertz_plastic_particle::ParticleType* particle2,
        std::chrono::duration<double>, const DEM::ParameterMap& parameters):
        psi0_(parameters.get_parameter<double>("psi0_")),
        psi0T_B_(parameters.get_parameter<double>("psi0T_B_")),
        kTp_(parameters.get_parameter<double>("kTp_")),
        kp_(parameters.get_parameter<double>("kp_")),
        R0_(parameters.get_parameter<double>("R0")),
        //Rb_(parameters.get_parameter<double>("Rb")),
        bt_(parameters.get_parameter<double>("bt")),
        br_(parameters.get_parameter<double>("br")),
        A(parameters.get_parameter<double>("A")),
        Ep_eff_(parameters.get_parameter<double>("Ep_eff")),
        kb_coeff(parameters.get_parameter<double>("kb_coeff")),
        binder_yield_stress_(parameters.get_parameter<double>("binder_yield_stress_")),
        yield_h_(parameters.get_parameter<double>("yield_h_")),
        h_(parameters.get_parameter<double>("h")),
        hmax_(parameters.get_parameter<double>("hmax")),
        h_plast_(parameters.get_parameter<double>("h_plast_")),
        mu_particle_(parameters.get_parameter<double>("mu_particle")),
        particle_yield_stress_(parameters.get_parameter<double>("particle_yield_stress_")), //comment here
        bonded_(parameters.get_parameter<bool>("bonded_")),
        adhesive_(parameters.get_parameter<bool>("adhesive_")),
        binder_contact_(parameters.get_parameter<bool>("binder_contact")),
        particle_contact_(parameters.get_parameter<bool>("particle_contact")),
        fractured_(parameters.get_parameter<bool>("fractured")),
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
                Bi.push_back(parameters.get_parameter<double>("A_" + std::to_string(i)));
                di_.push_back(parameters.get_parameter<double>("d_" + std::to_string(i)));
                ddi_.push_back(parameters.get_parameter<double>("dd_" + std::to_string(i)));
                dti_Scalar.push_back(parameters.get_parameter<double>("dti_Scalar"+ std::to_string(i)));
                ddti_Scalar.push_back(parameters.get_parameter<double>("ddti_Scalar"+ std::to_string(i)));
                dti_.push_back(parameters.get_vec3("dt_" + std::to_string(i)));
                ddti_.push_back(parameters.get_vec3("ddt_" + std::to_string(i)));
            }
        }

DEM::elastic_plastic_binder_hertz_plastic_particle::elastic_plastic_binder_hertz_plastic_particle(
        DEM::elastic_plastic_binder_hertz_plastic_particle::ParticleType* particle1,
        DEM::elastic_plastic_binder_hertz_plastic_particle::SurfaceType * surface1,
        std::chrono::duration<double>, const DEM::ParameterMap& parameters):
        psi0_(parameters.get_parameter<double>("psi0_")),
        psi0T_B_(parameters.get_parameter<double>("psi0T_B_")),
        kTp_(parameters.get_parameter<double>("kTp_")),
        kp_(parameters.get_parameter<double>("kp_")),
        R0_(parameters.get_parameter<double>("R0")),
        //Rb_(parameters.get_parameter<double>("Rb")),
        Ep_eff_(parameters.get_parameter<double>("Ep_eff")),
        kb_coeff(parameters.get_parameter<double>("kb_coeff")),
        binder_yield_stress_(parameters.get_parameter<double>("binder_yield_stress_")),
        yield_h_(parameters.get_parameter<double>("yield_h_")),
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
        particle_contact_(parameters.get_parameter<bool>("particle_contact")),
        fractured_(parameters.get_parameter<bool>("fractured")),
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
        Bi.push_back(parameters.get_parameter<double>("A_" + std::to_string(i)));
        di_.push_back(parameters.get_parameter<double>("d_" + std::to_string(i)));
        ddi_.push_back(parameters.get_parameter<double>("dd_" + std::to_string(i)));
        dti_Scalar.push_back(parameters.get_parameter<double>("dti_Scalar"+ std::to_string(i)));
        ddti_Scalar.push_back(parameters.get_parameter<double>("ddti_Scalar"+ std::to_string(i)));
        dti_.push_back(parameters.get_vec3("dt_" + std::to_string(i)));
        ddti_.push_back(parameters.get_vec3("ddt_" + std::to_string(i)));
    }
}

void DEM::elastic_plastic_binder_hertz_plastic_particle::update(double h, const DEM::Vec3& dt, const Vec3& drot, const DEM::Vec3& normal)
{
//    std::cout << "=======New iteration=======" << std::endl;
    rot_ += drot;
    F_ = update_normal_force(h);
//    std::cout <<uT_<<"," << dt<<","<< normal << std::endl;
    update_tangential_force(dt, normal);
//    std::cout << "h:"<< h << std::endl;
//      std::cout << "F_:"<< F_ << std::endl;
//     std::cout << "F_binder:"<< F_binder << std::endl;
//      std::cout << "F_particle:"<< F_particle << std::endl;
}

unsigned DEM::elastic_plastic_binder_hertz_plastic_particle::M;
const DEM::ElectrodeMaterial* DEM::elastic_plastic_binder_hertz_plastic_particle::material;

double  DEM::elastic_plastic_binder_hertz_plastic_particle::update_normal_force(double h)
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
        if (((h_ > -bt_) && !particle_contact_) || (bonded_ && !particle_contact_))
        {
            double viscoelastic_summation = 0.;
            for (unsigned i = 0; i != M; i++)
            {
                ddi_[i] = Ai[i] * ((bt_ + h_) - di_[i]) + Bi[i] * dh;
                viscoelastic_summation += alpha_i[i] * ddi_[i];
                di_[i] += ddi_[i];
            }
            F_binder += psi0_ * (dh - viscoelastic_summation);
//            if (F_binder <= 0 && !bonded_)
//            {
//                F_binder = 0;
//            }
            // Check yield criterion and remove the stress added if the material has yielded
            // Effective stress is set to uniaxial stress
            //double sigma_Mises_effective_binder = F_binder/(A*(1-v1));
            double sigma_binder = F_binder/A;
            if (sigma_binder > binder_yield_stress_)
            {
                F_binder = binder_yield_stress_*A;
            }
            else if (sigma_binder < -binder_yield_stress_)
            {
                F_binder = -binder_yield_stress_*A;
            }
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
//        if ((h_ > -bt_)  && adhesive() && !bonded_) //Change to if ((h_ > -bt_) && adhesive())?? kollar då om partiklarna är inom binder avstånd istället för på binderkraften
//        {
//            bonded_ = true;
//        }
        if ((F_binder > 0)  && adhesive() && !bonded_) //Change to if ((h_ > -bt_) && adhesive())?? kollar då om partiklarna är inom binder avstånd istället för på binderkraften
        {
            bonded_ = true;
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

// ===PARTICLE CONTACT MODEL============================================================================================
    if (h_ > 0)
    {
//        if (material->bond_breaking){
//            particle_contact_ = true;
//            bonded_ = false;
//        }
        if (h>=yield_h_ && h >= hmax_) // Perfectly plastic deformation of particles, new maximum overlap and contact radius
        {
            F_particle += 1.5 * kp_ * sqrt(yield_h_) * dh;
        }
        else if (h > h_plast_) //if overlap larger than elastic recover region
        {
            F_particle += 1.5*kp_*sqrt(h_)*dh;
            if(F_particle<0)
            {
                F_particle = 0;
                h_plast_ = h_;
            }
        }
        else
        {
            F_particle = 0;
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
    else
    {
        return std::max(F_particle, 0.) + std::max(F_binder, 0.);
    }
}

void  DEM::elastic_plastic_binder_hertz_plastic_particle::update_tangential_force(const DEM::Vec3& dt, const DEM::Vec3& normal)
{
    if (adhesive() && bonded_)
    {
//======NEW VISCOELASTIC MODEL IN TANGENITIAL DIRECTION, WORKING WITH SCALAR VALUES OF TANGENTIAL DISPLACEMENT==========
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
//======================================================================================================================

//======OLD TANGENTIAL FORCE UPDATE USING TANGENTIAL VECTOR==============================================================
//      FT_binder_ -= dot_product(FT_binder_,normal)*normal;
//      uT_ -= dot_product(uT_,normal)*normal;
//      uT_ += dt;
//      dFT_ = dt;
//      for (unsigned i = 0; i !=M; i++)
//      {
//          dti_[i] -= dot_product(dti_[i],normal)*normal;
//          ddti_[i] = Bi[i]*dt + Ai[i]*(uT_ - dti_[i]);
//          dFT_ -= alpha_i[i]*ddti_[i];
//          dti_[i] += ddti_[i];
//      }
//      FT_binder_ += psi0T_B_ * dFT_;
//    }
//======================================================================================================================

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

std::string DEM::elastic_plastic_binder_hertz_plastic_particle::restart_data() const {
        std::ostringstream ss;
        ss << named_print(dt_, "dt") << ", "
           << named_print(bt_, "bt") << ", "
           << named_print(br_, "br") << ", "
           << named_print(h_, "h") << ", "
           << named_print(hmax_, "hmax") << ", "
           << named_print(h_plast_, "h_plast_") << ", "
           << named_print(kTp_, "kTp_") << ", "
           << named_print(kp_, "kp_") << ", "
           << named_print(kb_coeff, "kb_coeff") << ", "
           << named_print(Ep_eff_, "Ep_eff") << ", "
           << named_print(binder_yield_stress_, "binder_yield_stress_") << ", "
           << named_print(particle_yield_stress_, "particle_yield_stress_") << ", "
           << named_print(yield_h_, "yield_h_") << ", "
           << named_print(R0_, "R0") << ", "
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
           << named_print(particle_contact_, "particle_contact") << ", "
           << named_print(fractured_, "fractured") << ", "
           << named_print(psi0_, "psi0_") << ", "
           << named_print(psi0T_B_, "psi0T_B_") << ", "
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

bool DEM::elastic_plastic_binder_hertz_plastic_particle::create_binder_contact(const ElectrodeMaterial* mat)
{
    std::random_device random_device;
    std::default_random_engine rand_engine { random_device() };
    std::uniform_real_distribution<double> distribution{0., 1.};
    double random_value = distribution(rand_engine);
    return random_value < mat->fraction_binder_contacts && mat->new_binder_contacts;
}

DEM::Vec3 DEM::elastic_plastic_binder_hertz_plastic_particle::get_rolling_resistance_torque() const {
    //if (binder_contact_) {
    //    return -Rb_*Rb_*0.01*kB_*rot_;
    //}
    //else {
    return {0, 0, 0};
    //    }
}

std::string DEM::elastic_plastic_binder_hertz_plastic_particle::get_output_string() const {
    std::stringstream ss;
    ss  << F_ << ", " << F_binder << ", " << F_particle << ", " << hmax_ << ", "
        << FT_.x() << ", "  << FT_.y() << ", " << FT_.z() << ", "
        << FT_binder_.x() << ", "  << FT_binder_.y() << ", " << FT_binder_.z() << ", "
        << FT_part_.x() << ", "  << FT_part_.y() << ", " << FT_part_.z() << ", "
        << binder_contact_ << ", " << bt_;
    return ss.str();
}

void DEM::elastic_plastic_binder_hertz_plastic_particle::set_increment(std::chrono::duration<double> dt) {
    dt_ = dt.count();
    Ai = {};
    Bi = {};
    for (unsigned i=0; i!=M; ++i) {
        Ai.push_back(1-exp((-dt_/tau_i[i])));
        Bi.push_back(tau_i[i]/dt_ *((dt_/tau_i[i])-Ai[i]));
    }
}

bool DEM::elastic_plastic_binder_hertz_plastic_particle::adhesive() const
{
    return adhesive_ && material->adhesive;
}