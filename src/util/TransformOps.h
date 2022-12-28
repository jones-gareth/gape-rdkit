//
// Created by jones on 12/23/2022.
//

#ifndef GAPE_TRANSFORMOPS_H
#define GAPE_TRANSFORMOPS_H

#include <Geometry/Transform3D.h>

namespace Gape {

    using namespace RDKit;

    /**
     * Determines the matrix transformation to perform a rotation (of angle
     * angle) about an arbitrary axis (defined by p1 and p2). The result is
     * returned in matrix.
     *
     * Ref: Newman & Sproull, Principles of Interactive computer Graphics,
     * McGraw Hill, 1981 pp 346-348.
     */
    void determineRotation(const RDGeom::Point3D &point1, const RDGeom::Point3D &point2, double angle,
                           RDGeom::Transform3D &transform);
} //namespace Gape
#endif //GAPE_TRANSFORMOPS_H
