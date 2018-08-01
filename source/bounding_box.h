//
// Created by erolsson on 2018-07-31.
//

#ifndef DEMSIM_BOUNDING_BOX_H
#define DEMSIM_BOUNDING_BOX_H

#include <algorithm>

#include "bounding_box_projection.h"
#include "surface.h"

namespace DEM {
    template <typename ForceModel, typename ParticleType>
    class BoundingBox {
    using BProjectionType = BoundingBoxProjection<ForceModel, ParticleType>;
    using SurfaceType = Surface<ForceModel, ParticleType>;
    public:
        explicit BoundingBox(ParticleType*);
        explicit BoundingBox(SurfaceType*);

        void update() {update_function(); };

        void set_stretch(double stretch) {stretch_ = stretch; }
        BProjectionType bx;
        BProjectionType ex;
        BProjectionType by;
        BProjectionType ey;
        BProjectionType bz;
        BProjectionType ez;

    private:
        ParticleType* particle_;
        SurfaceType* surface_;
        double stretch_{ 0. };

        // Function pointer to the update function, set in the construction of the bounding box
        // In doing so, if statements is avoided at each update
        void (BoundingBox<ForceModel, ParticleType>::*update_function)();
        void particle_update();
        void surface_update();
    };

    template<typename ForceModel, typename ParticleType>
    BoundingBox<ForceModel, ParticleType>::BoundingBox(ParticleType* p) :
        particle_(p),
        surface_(nullptr),
        update_function(&BoundingBox<ForceModel, ParticleType>::particle_update)
    {

    }

    template<typename ForceModel, typename ParticleType>
    BoundingBox<ForceModel, ParticleType>::BoundingBox(BoundingBox::SurfaceType* s) :
        particle_(nullptr),
        surface_(s),
        update_function(&BoundingBox<ForceModel, ParticleType>::surface_update)
    {

    }


    template<typename ForceModel, typename ParticleType>
    void BoundingBox<ForceModel, ParticleType>::particle_update()
    {
        Vec3 particle_position = particle_->get_position();
        double R = particle_->get_radius();
        bx.value = particle_position.x - R - stretch_;
        ex.value = particle_position.x + R + stretch_;

        by.value = particle_position.y - R - stretch_;
        ey.value = particle_position.y + R + stretch_;

        bz.value = particle_position.z - R - stretch_;
        ez.value = particle_position.z + R + stretch_;
    }

    template<typename ForceModel, typename ParticleType>
    void BoundingBox<ForceModel, ParticleType>::surface_update()
    {

    }

}
#endif //DEMSIM_BOUNDING_BOX_H
