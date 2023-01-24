//
// Created by gareth on 1/19/23.
//

#ifndef GAPE_PHARMFEATUREGEOMETRY_H
#define GAPE_PHARMFEATUREGEOMETRY_H

#include <Geometry/point.h>
#include <string>

namespace Gape {
    enum Geometry {
        Sphere, Vector, Arc, MultiVector
    };

    class PharmFeatureGeometry {
    protected:
        const Geometry geometry;
    public:
        virtual std::string summary() = 0;
        /**
         * Gets the number of points used by this feature.
         *
         * @return
         */
        virtual int getNumPoints() = 0;

        /**
         * Returns the coordinates of a point
         *
         * @param no
         * @return
         */
        virtual RDGeom::Point3D getPoint(int num) = 0;
    };

} // Gape

#endif //GAPE_PHARMFEATUREGEOMETRY_H
