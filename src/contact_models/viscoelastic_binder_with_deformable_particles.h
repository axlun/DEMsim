//
// Created by Axel on 2021-09-29.
//

#ifndef DEMSIM_VISCOELASTIC_BINDER_WITH_DEFORMABLE_PARTICLES_H
#define DEMSIM_VISCOELASTIC_BINDER_WITH_DEFORMABLE_PARTICLES_H

#include <chrono>
#include <vector>

#include "../particles/spherical_particle.h"
#include "../surfaces/surface_base.h"
#include "../surfaces/point_surface.h"
#include "../utilities/vec3.h"

namespace DEM {
    class ElectrodeMaterial;
    class ParameterMap;
    class Viscoelastic_binder_with_deformable_particles
    {
        using ParticleType = SphericalParticle<Viscoelastic_binder_with_deformable_particles>;
        using SurfaceType = Surface<Viscoelastic_binder_with_deformable_particles, ParticleType>;

    public:
        Viscoelastic_binder_with_deformable_particles(ParticleType *particle1, ParticleType *particle2,
                                                      std::chrono::duration<double> dt);

        [[nodiscard]] double get_overlap() const { return h_; }

        [[nodiscard]] double get_normal_force() const { return F_; }

        [[nodiscard]] bool active() const { return F_ != 0; }

        void update(double h, const Vec3& dt, const Vec3& drot, const Vec3& normal);


    private:
        const static ElectrodeMaterial* material;

        double kp_ = 0; //particle stiffness, inte noll men något bra
        //double kB_  =0;
        double psi0_ = 0; //initial stiffness,också nått bra här
        double dK_dhbn_ = 0;
        double Khbn_  = 0;
        double hb_ = 0;
        double hp_ = 0;
        double h_ = 0.;
        double F_ = 0;
        //double dF_ = 0;
        double bt_;
        double R0_;

        bool bonded_ = false;
        bool adhesive_ = true;
        bool binder_contact_ ;

        static unsigned M;
        double dt_;   // Time increment

        std::vector<double> tau_i{};
        std::vector<double> alpha_i{};
        std::vector<double> Ai{};
        std::vector<double> Bi{};
        std::vector<double> di_{};
        std::vector<double> ddi_{};
        std::vector <DEM::Vec3> ddti_{};
        std::vector <DEM::Vec3> dti_{};


        double update_normal_force(double h);
        static bool create_binder_contact(const ElectrodeMaterial* mat);
        bool adhesive() const;
    };
}

#endif //DEMSIM_VISCOELASTIC_BINDER_WITH_DEFORMABLE_PARTICLES_H
