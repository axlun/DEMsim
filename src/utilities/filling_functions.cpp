//
// Created by erolsson on 2018-09-10.
//

#include "filling_functions.h"
#include <vector>

#include <random>

std::vector<DEM::Vec3> DEM::random_fill_cylinder(double z0, double z1, double cylinder_radius,
        const std::vector<double>& radii, double distance_between_objects)
{
    std::vector<Vec3> particle_positions;
    std::random_device random_device;
    std::default_random_engine rand_engine(random_device());
    for (auto r : radii) {
        std::uniform_real_distribution<double> dist_r(-cylinder_radius+r + distance_between_objects,
                                                      cylinder_radius-r - distance_between_objects);
        std::uniform_real_distribution<double> dist_z(z0+r + distance_between_objects, z1-r - distance_between_objects);
        bool overlapping = true;
        Vec3 position {};
        while(overlapping) {
            position.x() = dist_r(rand_engine);
            position.y() = dist_r(rand_engine);
            position.z() = dist_z(rand_engine);

            //Check if a particle at the chosen position overlaps with an other
            if (position.x()*position.x()+position.y()*position.y() <
                   (cylinder_radius - r - distance_between_objects)*(cylinder_radius - r - distance_between_objects)) {
                overlapping = check_overlaps(position, r + distance_between_objects, particle_positions, radii);
            }
        }

        particle_positions.push_back(position);
    }
    return particle_positions;
}

bool DEM::check_overlaps(const DEM::Vec3& point, double radius, const std::vector<DEM::Vec3>& particle_positions,
                        const std::vector<double>& radii)
{
    for (unsigned i = 0; i != particle_positions.size(); ++i) {
        Vec3 position_i = particle_positions[i];
        if ((radii[i] + radius - (point-position_i).length()) > 0.) {
            return true;
        }
    }
    return false;
}

std::vector<DEM::Vec3> DEM::random_fill_box(double z0, double z1, double box_width,
                                           const std::vector<double>& radii, double bt)
{
    std::vector<Vec3> particle_positions;
    std::random_device random_device;
    std::default_random_engine rand_engine(random_device());
    for (auto r : radii) {
        std::uniform_real_distribution<double> dist_r(0.0+r + bt,
                                                      (box_width)-(r + bt ));

        std::uniform_real_distribution<double> dist_z(z0+(r+bt), z1-(r + bt));
        bool overlapping = true;
        Vec3 position {};
            while(overlapping) {
            position.x() = dist_r(rand_engine);
            position.y() = dist_r(rand_engine);
            position.z() = dist_z(rand_engine);
            //Check if a particle at the chosen position overlaps with an other
            overlapping = check_overlaps(position, r+bt, particle_positions, radii);
        }

        particle_positions.push_back(position);
    }
    return particle_positions;
}



