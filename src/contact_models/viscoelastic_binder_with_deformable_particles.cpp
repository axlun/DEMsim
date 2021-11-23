//
// Created by Axel on 2021-09-29.
//


#include "viscoelastic_binder_with_deformable_particles.h"

#include <cmath>
#include <iostream>
#include <random>
#include <vector>

#include "../materials/electrode_material.h"


DEM::Viscoelastic_binder_with_deformable_particles::Viscoelastic_binder_with_deformable_particles(DEM::Viscoelastic_binder_with_deformable_particles::ParticleType *particle1,
                                                        DEM::Viscoelastic_binder_with_deformable_particles::ParticleType *particle2,std::chrono::duration<double> dt)
{
    //extracting from material
    auto mat1 = dynamic_cast<const ElectrodeMaterial *>(particle1->get_material());
    auto mat2 = dynamic_cast<const ElectrodeMaterial *>(particle2->get_material());
    material = mat1;

    R0_= 1. / (1. / particle1->get_radius() + 1. / particle2->get_radius());

    double E1 = mat1->E;
    double v1 = mat1->nu;
    double E2 = mat1->E;
    double v2 = mat1->nu;
    double vp1 = mat1->nup;
    double vp2 = mat2->nup;
    double Ep2 = mat2->Ep;
    double Ep1 = mat1->Ep;
    double Ep_eff = 1./(((1-vp1*vp1)/Ep1)+((1-vp2*vp2)/Ep2));
    bt_ = mat1 -> bt; //Thickness of the binder link
    double br_ = mat1 -> binder_radius_fraction*2*R0_; //radius of the binder link, the factor 2 comes form using the effective radius and not the radius of the particle itself,
                                                        //this might be changed
    double A = DEM::pi*br_*br_; //area of the binder contact
    psi0_ = (1 - v1)/(1 + v1)/(1 - 2*v1)*E1*A; //The instantaneous value of the relaxation function for the binder material
    binder_contact_ = create_binder_contact(mat1);
    adhesive_ = true;
    M = mat1->M();
    tau_i = mat1->tau_i;
    alpha_i = mat1->alpha_i;

    dt_ = dt.count();// time increment
    kp_ = Ep1/(1-vp1*vp1)*br_;

    for (unsigned i=0; i!=M; ++i)
    {
        di_.push_back(0);
        ddi_.push_back(0);
        Ai.push_back((1-exp(-dt_/tau_i[i]))*(Khbn_-di_[i]));
        Bi.push_back((tau_i[i]/dt_) * dK_dhbn_ * ((dt_/tau_i[i])+exp(-dt_/tau_i[i])-1));
    }
}

void DEM::Viscoelastic_binder_with_deformable_particles::update(double h, const DEM::Vec3& dt, const Vec3& drot, const DEM::Vec3& normal)
{
    F_ = update_normal_force(h);
}

unsigned DEM::Viscoelastic_binder_with_deformable_particles::M;
const DEM::ElectrodeMaterial* DEM::Viscoelastic_binder_with_deformable_particles::material;

double DEM::Viscoelastic_binder_with_deformable_particles::update_normal_force(double h)
{
    std::cout << "**New iteration**" << std::endl;

    double dh = h - h_;
    h_ = h;
    std::cout << "h:" << h_ << std::endl;
    std::cout << "dh:" << dh << std::endl;


    if (binder_contact_)
    {
        if ((h > -bt_) || bonded_)
        {

            Khbn_ = hb_/(bt_-hb_);
            std::cout << "Khbn_:" << Khbn_ << std::endl;
            dK_dhbn_ = bt_/((bt_-hb_)*(bt_-hb_));
            std::cout << "dK_dhbn_:" << dK_dhbn_ << std::endl;
            for (unsigned i = 0; i != M; ++i)
            {
                Ai[i] = (1-exp(-dt_/tau_i[i]))*(Khbn_-di_[i]);
                //std::cout << "Khbn_:" << Khbn_<< std::endl;
                //std::cout << "di_[i]:" << di_[i]<< std::endl;
                std::cout << "Ai:" << Ai[i]<< std::endl;
                Bi[i] = (tau_i[i]/dt_) * dK_dhbn_ * ((dt_/tau_i[i])+exp(-dt_/tau_i[i])-1);
                std::cout << "Bi:" << Bi[i]<< std::endl;
            }
            double gamma = 0;
            double theta = 0;
            for (unsigned i = 0; i != M; ++i)
            {
                gamma -= alpha_i[i]*Bi[i];
                theta += alpha_i[i]*Ai[i];
            }
            gamma = (gamma + dK_dhbn_) * psi0_;
            std::cout << "Gamma:" << gamma << std::endl;
            theta = theta * psi0_;
            std::cout << "Theta:" << theta << std::endl;
            double dh_b = 0;
            double dh_p = 0;
            dh_b = (dh+theta/kp_)/(1+gamma/kp_);
            std::cout << "dh_b_:" << dh_b << std::endl;
            //dh_b = dh/(1+(gamma-theta)/kp_);
            dh_p = dh-dh_b;
            std::cout << "dh_p:" << dh_p << std::endl;
            hb_ = hb_+dh_b;
            std::cout << "hb_:" << hb_ << std::endl;
            hp_ = hp_+dh_p;
            std::cout << "hp_:" << hp_ << std::endl;

            std::cout << "h_:" << h_ << std::endl;


            for (unsigned i = 0; i != M; ++i)
            {
                ddi_[i] = Ai[i] + Bi[i]*dh_b;
                di_[i] += ddi_[i];
                std::cout << "di_[i]:" << di_[i]<< std::endl;
            }
            double dF_particle = kp_*dh_p;
            std::cout << "dF Particles:" << dF_particle << std::endl;

            double  visco_counter = 0;
            for (unsigned i = 0; i != M; ++i)
            {
                visco_counter += alpha_i[i]*(Ai[i]+Bi[i]*dh_b);
            }

            double dF_binder = psi0_*(dK_dhbn_*dh_b - visco_counter);
            std::cout << "dF Binder:" << dF_binder << std::endl;

            F_ = F_ + dF_particle;
            std::cout << "F_:" << F_ << std::endl;
        }
        else{
            F_ = 0;
            for (unsigned i = 0; i != M; ++i) {
                di_[i] = 0;
            }
        }
        if (F_ > 0 && adhesive()) {
            bonded_ = true;
        }
    }
    else{
        F_ = kp_*h_; //använd hertz istället
    }
    return F_;
}


bool DEM::Viscoelastic_binder_with_deformable_particles::create_binder_contact(const ElectrodeMaterial* mat)
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

bool DEM::Viscoelastic_binder_with_deformable_particles::adhesive() const {
    return adhesive_ && material->adhesive;
}