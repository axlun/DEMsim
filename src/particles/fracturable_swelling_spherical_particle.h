//
// Created by Axel on 2025-01-09.
//

#ifndef DEMSIM_FRACTURABLE_SWELLING_SPHERICAL_PARTICLE_H
#define DEMSIM_FRACTURABLE_SWELLING_SPHERICAL_PARTICLE_H

#include "particle_base.h"
#include "spherical_particle_base.h"
#include "../utilities/contact_vector.h"
#include "../materials/material_base.h"
#include "../engine/contact.h"
#include "../utilities/vec3.h"

namespace DEM
{

    template<typename ForceModel, typename FracturableSwellingSphericalParticle> class Contact;

    template<typename ForceModel>
    class FracturableSwellingSphericalParticle : public SphericalParticleBase<ForceModel>
    {
        using ContactType = Contact<ForceModel, FracturableSwellingSphericalParticle<ForceModel>>;
        using ContactPointerType = typename ContactMatrix<ContactType>::PointerType;

    public:
        // No assignment or plain copies of particles

        FracturableSwellingSphericalParticle(double radius, const Vec3& position, const Vec3& velocity,
                                 const MaterialBase* material, std::size_t object_id, std::size_t collision_id=0);
        FracturableSwellingSphericalParticle(const ParameterMap& parameters, MaterialBase* material);
        virtual ~FracturableSwellingSphericalParticle() = default;
        void sum_contact_forces()
        {
            SphericalParticleBase<ForceModel>::sum_contact_forces(contacts_);
        }

        void add_contact(ContactPointerType contact, std::size_t index, int direction)
        {
            SphericalParticleBase<ForceModel>::add_contact(contact, index, direction, contacts_);
        }

        void remove_contact(std::size_t index)
        {
            SphericalParticleBase<ForceModel>::remove_contact(index, contacts_);
        }
        [[nodiscard]] std::size_t get_numer_of_contacts() const {return contacts_.size();}
        ContactVector<std::pair<ContactPointerType, int>>& get_contacts() {return contacts_;}

        using SphericalParticleBase<ForceModel>::get_id;
        using SphericalParticleBase<ForceModel>::get_position;

        [[nodiscard]] double get_radius() const override { return swell_state_*radius_; }
        [[nodiscard]] double get_material_scaling() const { return material_scaling_; }
        void scale_material(const double new_material_scaling_this_inc) override;

        using SphericalParticleBase<ForceModel>::get_fracture;
        using SphericalParticleBase<ForceModel>::set_fracture;

        using SphericalParticleBase<ForceModel>::get_swell_rate;
        using SphericalParticleBase<ForceModel>::set_swell_rate;
        void swell(const double new_swelling_this_inc) override;

        using SphericalParticleBase<ForceModel>::kinetic_energy;

        [[nodiscard]] virtual std::string get_output_string() const override;
        [[nodiscard]] virtual std::string restart_data() const override;

    private:
        ContactVector<std::pair<ContactPointerType, int>> contacts_;
        using SphericalParticleBase<ForceModel>::material_;
        using SphericalParticleBase<ForceModel>::mass_;
        using SphericalParticleBase<ForceModel>::radius_;
        using SphericalParticleBase<ForceModel>::inertia_;
        double swell_state_{1.};
        using ParticleBase<ForceModel>::swell_rate_;
        using ParticleBase<ForceModel>::swelling_this_inc_;
        double material_scaling_{1.};
        using ParticleBase<ForceModel>::material_scale_rate_;
        using ParticleBase<ForceModel>::material_scaling_this_inc_;
//        bool fracture_{false};
    protected:
        using SphericalParticleBase<ForceModel>::id_;
        using SphericalParticleBase<ForceModel>::position_;
        using SphericalParticleBase<ForceModel>::rot_;
        using SphericalParticleBase<ForceModel>::f_;
    };

}

#include "fracturable_swelling_spherical_particle.tpp"

#endif //DEMSIM_SWELLING_SPHERICAL_PARTICLE_H
