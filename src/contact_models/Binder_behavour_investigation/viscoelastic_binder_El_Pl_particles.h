//
// Created by Axel on 2021-11-17.
//

#ifndef CATHODE_COMPACTION_CPP_VISCOELASTIC_BINDER_EL_PL_PARTICLES_H
#define CATHODE_COMPACTION_CPP_VISCOELASTIC_BINDER_EL_PL_PARTICLES_H
#include <chrono>
#include <vector>

#include "../../particles/spherical_particle.h"
#include "../../surfaces/surface_base.h"
#include "../../surfaces/point_surface.h"
#include "../../utilities/vec3.h"

namespace DEM{
    class ElectrodeMaterial;
    class ParameterMap;
    class viscoelastic_binder_El_Pl_particles{
        using ParticleType = SphericalParticle<viscoelastic_binder_El_Pl_particles>;
        using SurfaceType = Surface<viscoelastic_binder_El_Pl_particles, ParticleType>;

    public:
        viscoelastic_binder_El_Pl_particles(ParticleType* particle1, ParticleType* particle2, std::chrono::duration<double>dt);
        viscoelastic_binder_El_Pl_particles(ParticleType* particle1, SurfaceType* surface, std::chrono::duration<double>dt );

        viscoelastic_binder_El_Pl_particles(ParticleType*, ParticleType*, std::chrono::duration<double>, const ParameterMap& parameters);
        viscoelastic_binder_El_Pl_particles(ParticleType*, SurfaceType*, std::chrono::duration<double>, const ParameterMap& parameters);


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

        double kp_ = 0;
        double psi0_ = 0;
        double kb_coeff = 1;
        double kT_B_;
        double h_ = 0;
        double hmax_ = -1e99 ;
        double yield_h_ = 1e99;

        double bt_ = 0;
        double br_ = 0;
        double R0_ = 0;
        bool bonded_ = false;
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
        std::vector <DEM::Vec3> ddti_{};
        std::vector <DEM::Vec3> dti_{};

        double dF_ = 0.;
        double F_ = 0;
        double F_binder = 0;
        double F_particle = 0;
        double mu_particle_;
        Vec3 dFT_{Vec3(0., 0., 0.)};
        Vec3 FT_{Vec3(0., 0., 0.)};
        Vec3 FT_binder_ {Vec3(0., 0., 0.)};
        Vec3 FT_part_ {Vec3(0., 0., 0.)};
        Vec3 uT_{ Vec3(0., 0., 0.) };

        Vec3 rot_ {Vec3(0., 0., 0.)};


        double update_normal_force(double h);
        static bool create_binder_contact(const ElectrodeMaterial* mat);
        bool adhesive() const;
    };
}





#endif //CATHODE_COMPACTION_CPP_VISCOELASTIC_BINDER_EL_PL_PARTICLES_H
