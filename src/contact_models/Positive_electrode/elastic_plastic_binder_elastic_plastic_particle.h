//
// Created by Axel on 2022-10-16.
//

#ifndef ELASTIC_PLASTIC_BINDER_ELASTIC_PLASTIC_PARTICLE_H
#define ELASTIC_PLASTIC_BINDER_ELASTIC_PLASTIC_PARTICLE_H
#include <chrono>
#include <vector>

#include "../../particles/spherical_particle.h"
#include "../../surfaces/surface_base.h"
#include "../../surfaces/point_surface.h"
#include "../../utilities/vec3.h"

namespace DEM{
    class ElectrodeMaterial;
    class ParameterMap;
    class elastic_plastic_binder_elastic_plastic_particle{
        using ParticleType = SphericalParticle<elastic_plastic_binder_elastic_plastic_particle>;
        using SurfaceType = Surface<elastic_plastic_binder_elastic_plastic_particle, ParticleType>;

    public:
        elastic_plastic_binder_elastic_plastic_particle(ParticleType* particle1, ParticleType* particle2, std::chrono::duration<double>dt);
        elastic_plastic_binder_elastic_plastic_particle(ParticleType* particle1, SurfaceType* surface, std::chrono::duration<double>dt );

        elastic_plastic_binder_elastic_plastic_particle(ParticleType*, ParticleType*, std::chrono::duration<double>, const ParameterMap& parameters);
        elastic_plastic_binder_elastic_plastic_particle(ParticleType*, SurfaceType*, std::chrono::duration<double>, const ParameterMap& parameters);


        [[nodiscard]] double get_overlap() const {return h_;};
        [[nodiscard]] double get_normal_force() const {return F_;};
        [[nodiscard]] const Vec3& get_tangential_force() const { return FT_; }
        [[nodiscard]] Vec3 get_rolling_resistance_torque() const;
        [[nodiscard]] double active() const {return F_ !=0; };
        [[nodiscard]] std::string restart_data() const;
        [[nodiscard]] std::string get_output_string() const;

        void update(double h, const Vec3& dt, const Vec3& drot, const Vec3& normal);
        void set_increment(std::chrono::duration<double>);

    private:
        const static ElectrodeMaterial* material;

        double F_1_ = 0.0;
        double alpha_1_ = 0.0;
        double F_2_ = 0.0;
        double alpha_2_ = 0.0;
        double a_1_ = 0.0;
        double beta_1_ = 0.0;
        double a_2_ = 0.0;
        double beta_2_ = 0.0;

        double kp_ = 0.;
        double kTp_ = 0.;
        double psi0_ = 0;
        double psi0T_B_ = 0.;
        double h_ = 0.0;
        double hmax_ = -1e99;
        double h_plast_ =0.0;
        double a_0_ = -1e99;
        double F_0_ = -1e99; //Maximum particle contact force
        double particle_yield_stress_ = 1e99; //high stress if yield is not used
        double binder_yield_stress_ = 1e99; //high stress if yield is not used
        double A = 0.0;
        double v1 = 0.0;
        double bt_ = 0;
        double br_ = 0;
        double R0_ = 0;
        double Ep_eff_ = 0;
        bool bonded_ = false;
        bool particle_contact_ = false;
        bool adhesive_ = false;
        bool binder_contact_ = false;
        bool fractured_ = true;

        static unsigned M;
        double dt_;

        std::vector<double> tau_i{};
        std::vector<double> alpha_i{};
        std::vector<double> Ai{};
        std::vector<double> Bi{};
        std::vector<double> di_{};
        std::vector<double> ddi_{};
        std::vector<double> dti_Scalar{};
        std::vector<double> ddti_Scalar{};
        std::vector <DEM::Vec3> ddti_{};
        std::vector <DEM::Vec3> dti_{};

        double dF_ = 0.;
        double F_ = 0;
        double F_binder = 0;
        double FT_binder_Scalar_ = 0;
        double F_particle = 0;
        double mu_particle_ = 0.;
        Vec3 dFT_{Vec3(0., 0., 0.)};
        Vec3 FT_{Vec3(0., 0., 0.)};
        Vec3 FT_binder_ {Vec3(0., 0., 0.)};
        Vec3 FT_part_ {Vec3(0., 0., 0.)};
        Vec3 uT_{ Vec3(0., 0., 0.) };
        Vec3 rot_ {Vec3(0., 0., 0.)};

        double update_normal_force(double h);
        void update_tangential_force(const Vec3& dt, const Vec3& normal);
        static bool create_binder_contact(const ElectrodeMaterial* mat);
        bool adhesive() const;
    };
}





#endif //ELASTIC_PLASTIC_BINDER_ELASTIC_PLASTIC_PARTICLE_H