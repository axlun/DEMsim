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
                                                        DEM::Viscoelastic::ParticleType* particle2,std::chrono::duration<double> dt)
{
    //extracting from material
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
    bt_ = mat1 -> bt; //Thickness of the binder link
}

double DEM::Viscoelastic_binder_with_deformable_particles::update_normal_force(double h)
{
    double dh = h - h_;
    h_ = h;

    if (binder_contact_) {
        if ((h > -bt_) || bonded_)
        {
            double gamma = 0;
            double theta = 0;
            for (unsigned i = 0; i != M; ++i)
            {
                gamma -= alpha_i[i]*Bi[i];
                theta += alpha_i[i]*Ai[i];
            }
            gamma = (gamma + dK_dhbn_) * psi0_;
            theta = theta * psi0_;
            double dh_b = 0;
            double dh_p = 0;
            dh_b = dh/(1+(gamma-theta)/kp_);
            dh_p = dh-dh_b;
            hb_ = hb_+dh_b;
            hp_ = hp_+dh_p;

            //double dF = psi0_*dK_dhbn_*dh_b;
            for (unsigned i = 0; i != M; ++i)
            {
                ddi_[i] = Ai[i] + Bi[i]*dh_b;
                di_[i] += ddi_[i];
                Ai.push_back((1-exp(-dt_/tau_i))*(Khbn_-di_)); // update [i] not pushback
                Bi.push_back((tau_i/dt_) * dK_dhnb * ((dt_/tau_i)+exp(-dti_/tau_i)-1));
            }
            dK_dhbn_ = bt_/((bt_-hb_)*(bt_-hb_));

            double dF = kp_*dh_p;
            F_ = F_ + dF;
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
        F_ = kp_*h_; //hertz
    }

}


bool DEM::Viscoelastic::create_binder_contact(const ElectrodeMaterial* mat) {
    std::random_device random_device;
    std::default_random_engine rand_engine { random_device() };
    std::uniform_real_distribution<double> distribution{0., 1.};
    double random_value = distribution(rand_engine);
    if (random_value < mat->fraction_binder_contacts){
        return true;
    }
    return false;
}