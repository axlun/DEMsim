//
// Created by  Axel on 2024-01-24
//

#include "swelling_spherical_particle.h"


template< typename ForceModel>
DEM::SwellingSphericalParticle<ForceModel>::SwellingSphericalParticle(
        double radius, const Vec3& position,
        const Vec3& velocity,
        const MaterialBase* material,
        std::size_t object_id,
        std::size_t collision_id)
        :SphericalParticleBase<ForceModel>(radius, position, velocity, material, object_id, collision_id)
{
    //Empty constructor
}

template <typename ForceModel>
DEM::SwellingSphericalParticle<ForceModel>::SwellingSphericalParticle(
        const DEM::ParameterMap& parameters, DEM::MaterialBase *material)
        :SphericalParticleBase<ForceModel>(parameters, material),
        swell_state_(parameters.get_parameter<double>("swell_state"))
        /*,
        *swell_rate_(parameters.get_parameter<double>("swell_rate"))
        */
{
    //Empty constructor
}
template <typename ForceModel>
void DEM::SwellingSphericalParticle<ForceModel>::swell(const double new_swelling_this_inc)
{
    swelling_this_inc_ = new_swelling_this_inc;
    swell_state_ += swelling_this_inc_;
    inertia_ = 0.4 * mass_ * swell_state_ * radius_ * swell_state_ * radius_;
}


template <typename ForceModel>
std::string DEM::SwellingSphericalParticle<ForceModel>::get_output_string() const
{
    std::ostringstream ss;
    ss << id_ << ", " << position_.x() << ", " << position_.y() << ", " << position_.z() << ", "
       << rot_.x() << ", " << rot_.y() << ", " << rot_.z() << ", " << radius_ * swell_state_ << ", "
       <<  kinetic_energy() << ", " << material_->id << ", "
       << f_.x() << ", " << f_.y() << ", " << f_.z();
    return ss.str();
}

template <typename ForceModel>
std::string DEM::SwellingSphericalParticle<ForceModel>::restart_data() const
{
    using DEM::named_print;
    std::ostringstream ss;
    ss << SphericalParticleBase<ForceModel>::restart_data() << ", "
       << named_print(swell_state_, "swell_state");
       //other variables needed?
    return ss.str();
}
/*
* template<typename ForceModel>
* void DEM::SwellingSphericalParticle<ForceModel>::move(const DEM::Vec3& new_disp_this_inc)
* {
*     swell_state_ += swell_rate_;//dt; //How will this be fixed?
*     SphericalPartileBase<ForceModel>::move(new_disp_this_inc);
* }
*/


