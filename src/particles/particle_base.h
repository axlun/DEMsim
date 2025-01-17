//
// Created by erolsson on 2018-07-29.
//

#ifndef DEMSIM_PARTICLE_BASE_NEW_H
#define DEMSIM_PARTICLE_BASE_NEW_H

#include <memory>
#include <string>
#include <vector>

#include "../materials/material_base.h"
#include "../utilities/vec3.h"

namespace DEM {
    class ParameterMap;

    template<typename ForceModel>
    class ParticleBase {
    public:
        ParticleBase& operator=(const ParticleBase&) = delete;
        ParticleBase(double, const Vec3&, const Vec3&, const MaterialBase*, unsigned );
        ParticleBase(const ParameterMap& parameters, MaterialBase* material);
        virtual ~ParticleBase() = default;
        [[nodiscard]] unsigned get_id() const { return id_; }
        [[nodiscard]] const MaterialBase* get_material() const { return material_; }
        [[nodiscard]] const Vec3& get_force() const { return f_; }
        [[nodiscard]] const Vec3& get_torque() const { return torque_; }

        [[nodiscard]] double get_mass() const { return mass_; }

        virtual void swell(const double dt){};
        virtual void scale_material(const double dt){};

        [[nodiscard]] const Vec3& get_position() const { return position_; }
        [[nodiscard]] const Vec3& get_velocity() const { return velocity_; }
        [[nodiscard]] const Vec3& get_acceleration() const { return acceleration_; }
        void set_velocity(const Vec3& new_velocity) { velocity_ = new_velocity; }
        void set_acceleration(const Vec3& new_acceleration) { acceleration_ = new_acceleration; }

        [[nodiscard]] const Vec3& get_rotation() const { return rot_; }
        [[nodiscard]] const Vec3& get_angular_velocity() const { return ang_vel_; }
        [[nodiscard]] const Vec3& get_angular_acceleration() const {return angular_acceleration_; }
        void set_angular_velocity(Vec3 new_ang_vel) { ang_vel_ = new_ang_vel; }
        void set_angular_acceleration(Vec3 new_ang_acc) { angular_acceleration_ = new_ang_acc; }

        [[nodiscard]] const Vec3& get_displacement_this_increment() const { return displacement_this_inc_; }
        [[nodiscard]] const Vec3& get_rotation_this_increment() const { return rot_this_inc_; }

        [[nodiscard]] double get_swelling_this_increment() const {return swelling_this_inc_;}
        [[nodiscard]] double get_swell_rate() const {return swell_rate_;}
        void set_swell_rate(const double& swell_rate) {swell_rate_ = swell_rate;}

        [[nodiscard]] double get_material_scaling_this_increment() const {return material_scaling_this_inc_;}
        [[nodiscard]] double get_material_scale_rate() const {return material_scale_rate_;}
        void set_material_scale_rate(const double& material_scale_rate) {material_scale_rate_ = material_scale_rate;}

        [[nodiscard]] bool get_fracture() const {return fracture_;}
        void set_fracture() {fracture_ = true;}

        [[nodiscard]] virtual std::string restart_data() const;
        void reset_contact_forces() {f_.set_zero(); torque_.set_zero(); }
    protected:
        const unsigned id_;
        double mass_;
        const MaterialBase* material_;

        Vec3 position_;
        Vec3 velocity_;
        Vec3 acceleration_ {0, 0, 0};
        Vec3 rot_{ Vec3(0., 0., 0.) };
        Vec3 ang_vel_{ Vec3(0., 0., 0.) };
        Vec3 angular_acceleration_ {0, 0, 0};

        // Forces
        Vec3 f_{ Vec3(0., 0., 0.) };
        Vec3 torque_{ Vec3(0., 0., 0.) };

        double swelling_this_inc_{0};
        double swell_rate_{0};

        double material_scaling_this_inc_{0};
        double material_scale_rate_{0};

        bool fracture_ {false};

        // Needed for sticking friction model
        Vec3 displacement_this_inc_{ Vec3(0., 0., 0.)};
        Vec3 rot_this_inc_{ Vec3(0., 0., 0.)};

        std::vector<Vec3> contact_forces_ { std::vector<Vec3>() };
        unsigned number_of_contacts_ { 0 };
    };
}

#include "particle_base.tpp"

#endif //DEMSIM_PARTICLE_BASE_NEW_H