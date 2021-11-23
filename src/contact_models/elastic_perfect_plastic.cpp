//
// Created by Axel on 2021-11-15.
//

#include "elastic_perfect_plastic.h"

#include <cmath>
#include <iostream>
#include <random>
#include <vector>

#include "../materials/electrode_material.h"

DEM::Elastic_perfectly_plastic_binder::Elastic_perfectly_plastic_binder(DEM::Elastic_perfectly_plastic_binder::ParticleType *particle1,
                                                                        DEM::Elastic_perfectly_plastic_binder::ParticleType *particle2,
                                                                        std::chrono::duration<double> dt)
{
  //Extracting material
  auto mat1 = dynamic_cast<const ElectrodeMaterial *>(particle1->get_material());
  auto mat2 = dynamic_cast<const ElectrodeMaterial *>(particle2->get_material());
  material = mat1;

  R0_ = 1./(1./particle1->get_radius()+1./particle2->get_radius()); //Effective radius of the two particles
  double E1 = mat1->E;
  double v1 = mat1->nu;
  double E2 = mat2->E;
  double v2 = mat2->E;
  double Ep1 = mat1->Ep;
  double Ep2 = mat2->Ep;
  double vp1 = mat1->nup;
  double vp2 = mat2->nup;
  double sigma_b_yield = mat1->Syb;
  double Ep_eff = 1./(((1-vp1*vp1)/Ep1)+((1-vp2*vp2)/Ep2));
  bt_ = mat1 -> bt; //total thickness of binder link
  double br_ = mat1 -> binder_radius_fraction*particle1->get_radius(); //Binder radius following the particle radius
  double A = DEM::pi*br_*br_; //area of binder contact
  kb_ = (1 - v1)/(1 + v1)/(1 - 2*v1)*E1*A/bt_; //stiffness of binder from Peter and Pelles analytical model, small deformation
  kp_ = (4./3)*Ep_eff* sqrt(R0_); //Particle stiffness following Hertz contact for two particles
  h_binder_yield_ = sigma_b_yield*A/kb_-bt_; //Overlap which give failure in binder
  binder_contact_ = create_binder_contact(mat1);
}

void DEM::Elastic_perfectly_plastic_binder::update(double h, const DEM::Vec3& dt, const Vec3& drot, const DEM::Vec3& normal)
{
    F_ = update_normal_force(h);
}

unsigned DEM::Elastic_perfectly_plastic_binder::M;
const DEM::ElectrodeMaterial* DEM::Elastic_perfectly_plastic_binder::material;


double DEM::Elastic_perfectly_plastic_binder::update_normal_force(double h)
{
    //std::cout<<"**New iteration**" <<std::endl;

    double dh = h-h_;
    h_= h;

    if (binder_contact_ && h_<0.) //Bindercontact and the particles are not in contact
    {
        if(h>-bt_ || bonded_) //The particles are closer to another then the thickness of the binder or the particles are bonded.
        {

            if(h > h_binder_yield_) //Check if the binder has collapsed, maybe add condition which decreases binder thickness when plastic deformation appears
            {
                F_ = kb_*(bt_+h_binder_yield_);
            }
            else //if not use linear model
            {
                F_ =kb_*(bt_+h_);
            }
        }
        else //if no binder contact the force is zero
        {
            F_ = 0;
        }
        if (F_>0 && adhesive()) //when there is compression of the binder the particles become bonded
        {
            bonded_ = true;
            bonded_history_ = true;
        }
    }
    else //If there is no binder contact or the particles are in contact
    {
        F_ = kp_*pow(h,1.5);
        if (F_< kb_*(bt_+h_binder_yield_) && dh >=0 && bonded_history_) //if the particle contact force is lower than the binder collapse force
            // and the particles are going towards each other and the particles have been bonded before
        {
            F_ =kb_*(bt_+h_binder_yield_);
        }

        bonded_ = false;
        binder_contact_ = false;

        if(h_<0)
        {
            F_ = 0;
        }
    }
//    std::cout << "F_:" << F_ << std::endl;
    return F_;

}

bool DEM::Elastic_perfectly_plastic_binder::create_binder_contact(const ElectrodeMaterial* mat) {
    std::random_device random_device;
    std::default_random_engine rand_engine { random_device() };
    std::uniform_real_distribution<double> distribution{0., 1.};
    double random_value = distribution(rand_engine);
    if (random_value < mat->fraction_binder_contacts){
        return true;
    }
    return false;
}

bool DEM::Elastic_perfectly_plastic_binder::adhesive() const {
    return adhesive_ && material->adhesive;
}