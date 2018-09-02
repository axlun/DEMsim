//
// Created by erolsson on 2018-07-31.
//

#ifndef DEMSIM_COLLISION_DETECTOR_H
#define DEMSIM_COLLISION_DETECTOR_H

#include <vector>
#include <omp.h>

#include "bounding_box.h"
#include "bounding_box_projection.h"
#include "contact_matrix.h"
#include "contact_vector.h"
#include "cylinder.h"
#include "point_surface.h"

namespace DEM {
    template <typename ForceModel, typename ParticleType> class Contact;
    template <typename ForceModel, typename ParticleType>
    class CollisionDetector {
        using ContactType = Contact<ForceModel, ParticleType>;
    public:

        using BoundingBoxType = BoundingBox<ForceModel, ParticleType>;
        // using CollisionPair = std::pair<const BoundingBoxType*, const BoundingBoxType*>;
        using BoundingBoxProjectionType = BoundingBoxProjection<ForceModel, ParticleType>;
        using SurfaceType = Surface<ForceModel, ParticleType>;
        using CylinderType = Cylinder<ForceModel, ParticleType>;

        class CollisionPair {
        public:
            ParticleType* particle1;
            ParticleType* particle2;
            SurfaceType* surface;

            CollisionPair(ParticleType* p1, ParticleType* p2) :
                particle1(p1), particle2(p2), surface(nullptr)
            {
                // Empty constructor
            }

            CollisionPair(ParticleType* p1, SurfaceType* surf) :
                    particle1(p1), particle2(nullptr), surface(surf)
            {
                // Empty constructor
            }

            std::pair<std::size_t, size_t> get_id_pair() const {
                if (surface == nullptr)
                    return std::make_pair(particle1->get_id(), particle2->get_id());
                else
                    return std::make_pair(particle1->get_id(), surface->get_id());
            }
        };

        CollisionDetector(const std::vector<ParticleType*>& particles,
                          const std::vector<SurfaceType*>& surfaces,
                          const ContactMatrix<ContactType>& contacts);

        void setup();
        void do_check();  //Not const due to re-ordering of the proj vectors
        std::vector<CollisionPair> contacts_to_create() const { return contacts_to_create_.get_objects(); }
        std::vector<CollisionPair> contacts_to_destroy() const { return contacts_to_destroy_;}

    private:
        std::vector<BoundingBox<ForceModel, ParticleType> > bounding_boxes_{};

        std::vector<BoundingBoxProjectionType*> xproj_{};
        std::vector<BoundingBoxProjectionType*> yproj_{};
        std::vector<BoundingBoxProjectionType*> zproj_{};

        std::size_t n_ = 0;

        const std::vector<ParticleType*>& particles_;
        const std::vector<SurfaceType*>& surfaces_;
        const ContactMatrix<ContactType>& contacts_;

        ContactVector<CollisionPair, std::pair<std::size_t, std::size_t>> contacts_to_create_{};
        std::vector<CollisionPair> contacts_to_destroy_ {};

        //BoolMatrix activeCollisions;
        void update_bounding_boxes();

        void check_bounding_box_vector(std::vector<BoundingBoxProjectionType*>& vector);
        bool check_other_axes(const BoundingBoxProjectionType* b1, const BoundingBoxProjectionType* b2) const;

    };

    template<typename ForceModel, typename ParticleType>
    void CollisionDetector<ForceModel, ParticleType>::do_check()
    {
        contacts_to_create_.clear();
        contacts_to_destroy_.clear();
        update_bounding_boxes();
        check_bounding_box_vector(xproj_);
        check_bounding_box_vector(yproj_);
        check_bounding_box_vector(zproj_);
    }

    template<typename ForceModel, typename ParticleType>
    CollisionDetector<ForceModel, ParticleType>::CollisionDetector(const std::vector<ParticleType*>& particles,
                                                                   const std::vector<SurfaceType*>& surfaces,
                                                                   const ContactMatrix<ContactType>& contacts) :
        particles_(particles),
        surfaces_(surfaces),
        contacts_(contacts)
    {
        //Empty constructor
    }

    template<typename ForceModel, typename ParticleType>
    void CollisionDetector<ForceModel, ParticleType>::setup()
    {
        std::size_t counter = 0;
        // bounding_boxes_.reserve(particles_.size() + surfaces_.size());
        for(const auto& p: particles_){
            bounding_boxes_.emplace_back(p, counter);
            ++counter;
        }

        for(const auto& s: surfaces_){
            auto cylinder_ptr = dynamic_cast<CylinderType*>(s);
            if (cylinder_ptr != nullptr) {
                // Inspect the cylinder bounding box to figure out if it is inward or not
                const auto& cyl_bounding_box = cylinder_ptr->get_bounding_box_values();
                if (std::abs(cyl_bounding_box[1]-cyl_bounding_box[0]) < 2*cylinder_ptr->get_radius())
                    bounding_boxes_.emplace_back(cylinder_ptr, counter, true);
                else
                    bounding_boxes_.emplace_back(s, counter);
            }
            else {
                bounding_boxes_.emplace_back(s, counter);
            }
            ++counter;
        }

        for(auto& bounding_box: bounding_boxes_){
            xproj_.push_back(&bounding_box.bx);
            xproj_.push_back(&bounding_box.ex);

            yproj_.push_back(&bounding_box.by);
            yproj_.push_back(&bounding_box.ey);

            zproj_.push_back(&bounding_box.bz);
            zproj_.push_back(&bounding_box.ez);
        }
        n_ = xproj_.size();
    }

    template<typename ForceModel, typename ParticleType>
    void CollisionDetector<ForceModel, ParticleType>::update_bounding_boxes()
    {
        #pragma omp parallel for
        for(std::size_t i = 0; i < bounding_boxes_.size(); ++i){
            bounding_boxes_[i].update();
        }
    }

    template<typename ForceModel, typename ParticleType>
    void CollisionDetector<ForceModel, ParticleType>::check_bounding_box_vector(
            std::vector<CollisionDetector::BoundingBoxProjectionType*>& vector)
    {
        for (unsigned i = 0; i != n_; ++i){
             unsigned j = i;
             while (j != 0 && vector[j-1]->get_value() > vector[j]->get_value()){

                 BoundingBoxProjectionType* bbm = vector[j];
                 BoundingBoxProjectionType* bbn = vector[j-1];

                 //depending on de beginnings and endings of the swapping
                 // remove or add contact
                 char c1 = bbm->get_position_char();
                 char c2 = bbn->get_position_char();

                 if (c1 == 'e' && c2 == 'b') {
                     auto id_pair = std::make_pair(bbm->get_id(), bbn->get_id());
                     if (!contacts_to_create_.erase(id_pair)) {
                         if (contacts_.exist(id_pair.first, id_pair.second)) {
                             // There is actually a contact to destroy
                             if (bbm->get_particle() != nullptr && bbn->get_particle() !=nullptr)
                                 contacts_to_destroy_.push_back(CollisionPair(bbm->get_particle() ,bbn->get_particle()));
                             else if (bbm->get_particle() != nullptr && bbn->get_surface() != nullptr)
                                 contacts_to_destroy_.push_back(CollisionPair(bbm->get_particle() ,bbn->get_surface()));
                             else if (bbn->get_particle() != nullptr && bbm->get_surface() != nullptr)
                                 contacts_to_destroy_.push_back(CollisionPair(bbn->get_particle() ,bbm->get_surface()));
                         }
                     }
                 }
                 else if (c1 == 'b' && c2 == 'e') {
                     if (check_other_axes(bbm, bbn)) {
                         if (bbm->get_particle() != nullptr && bbn->get_particle() !=nullptr)
                             contacts_to_create_.insert(std::make_pair(bbm->get_id(), bbn->get_id()),
                                     CollisionPair(bbm->get_particle(), bbn->get_particle()));
                         else if (bbm->get_particle() != nullptr && bbn->get_surface() != nullptr)
                             contacts_to_create_.insert(std::make_pair(bbm->get_id(), bbn->get_id()),
                                     CollisionPair(bbm->get_particle(), bbn->get_surface()));
                         else if (bbn->get_particle() != nullptr && bbm->get_surface() != nullptr)
                             contacts_to_create_.insert(std::make_pair(bbn->get_id(), bbm->get_id()),
                                     CollisionPair(bbn->get_particle(), bbm->get_surface()));
                     }
                 }

                 std::swap(vector[j], vector[j-1]);

                 bbn->increase_index();
                 bbm->decrease_index();
                 --j;
             }
        }
    }


    template<typename ForceModel, typename ParticleType>
    bool
    CollisionDetector<ForceModel, ParticleType>::check_other_axes(const CollisionDetector::BoundingBoxProjectionType* b1,
            const CollisionDetector::BoundingBoxProjectionType* b2) const
    {
        auto idx1 = b1->get_indices_on_other_axes();
        auto idx2 = b2->get_indices_on_other_axes();

        // checking the first of the axes
        if ( (*idx1[0] < *idx2[0] && *idx2[0] < *idx1[1]) || (*idx2[0] < *idx1[0] && *idx1[0] < *idx2[1]) ) {
            // Check the second axis
            if ( (*idx1[2] < *idx2[2] && *idx2[2] < *idx1[3]) || (*idx2[2] < *idx1[2] && *idx1[2] < *idx2[3]) ) {
                return true;
            }
            return false;
        }
        return false;
    }

}

#endif //DEMSIM_COLLISION_DETECTOR_H
