//
// Created by Axel on 2021-11-15.
//

#ifndef DEM_ELASTIC_PERFECT_PLASTIC_H
#define DEM_ELASTIC_PERFECT_PLASTIC_H

#include <chrono>
#include <vector>

#include "../particles/spherical_particle.h"
#include "../surfaces/surface_base.h"
#include "../surfaces/point_surface.h"
#include "../utilities/vec3.h"

namespace DEM{
    class ElectrodeMaterial;
    class ParameterMap;
    class Elastic_perfectly_plastic_binder{
            using ParticleType = SphericalParticle<Elastic_perfectly_plastic_binder>;
            using SurfaceType = Surface<Elastic_perfectly_plastic_binder, ParticleType>;
    public:
        Elastic_perfectly_plastic_binder(ParticleType *particle1, ParticleType *particle2, std::chrono::duration<double> dt);

        [[nodiscard]] double get_overlap() const {return h_;}

        [[nodiscard]] double get_normal_force() const {return F_; }

        [[nodiscard]] bool active() const {return F_ != 0; }

        void update(double h, const Vec3& dt, const Vec3& drot, const Vec3& normal);

    private:
        const static ElectrodeMaterial* material;
        double kp_=0.; //particle stiffness
        double kb_ = 0.; //binder stiffness
        double hb_ = 0;
        double hp_ = 0;
        double h_ = 0.;
        double F_ = 0;
        double bt_ = 0;
        double R0_ = 0;
        double h_binder_yield_ = 0;

        bool bonded_ = false;
        bool bonded_history_ = false;
        bool adhesive_ = true;
        bool binder_contact_ = false;

        static unsigned M;
        double dt_; //time increment

        double update_normal_force(double h);
        static bool create_binder_contact(const ElectrodeMaterial* mat);
        bool adhesive() const;
    };

}

#endif //DEM_ELASTIC_PERFECT_PLASTIC_H
