/*
 * ShapeMetrics.h
 *
 *  Created on: Apr 12, 2016
 *      Author: gjones
 *
 * See https://en.wikipedia.org/wiki/Gyration_tensor for information
 * about the gyration tensor and accosiated shape descriptors
 */

#ifndef SRC_UTIL_SHAPEMETRICS_H_
#define SRC_UTIL_SHAPEMETRICS_H_

#include <Eigen/Dense>
#include "CoordOps.h"

namespace Gape {

using namespace std;
using namespace Eigen;

/**
 * A class for determining shape descriptors from the gyration tensor.
 * The gyration tensor describes the second moments of postions for atoms in a molecule.
 * Diagonalization yields three orientation independent descriptors/principal moments
 *
 * See https://en.wikipedia.org/wiki/Gyration_tensor
 */
class ShapeMetrics {
public:
    ShapeMetrics(const CoordMatrix & c) :
            coords(c) {
        calculateGyrationTensorShapeDescriptors();
    }

    virtual ~ShapeMetrics() {
    }

    ShapeMetrics(const ShapeMetrics & rhs) = delete;
    ShapeMetrics & operator =(const ShapeMetrics & rhs) = delete;
    ShapeMetrics(ShapeMetrics && rhs) = delete;
    ShapeMetrics & operator =(ShapeMetrics && rhs) = delete;

    static const Matrix3d determineGyrationTensor(const CoordMatrix & coords);

    // getters to retieve shape descriptors

    const double getAcylindricity() const {
        return acylindricity;
    }

    const double getAnisotrophy() const {
        return anisotrophy;
    }

    const double getAsphericity() const {
        return asphericity;
    }

    const CoordMatrix & getCoords() const {
        return coords;
    }

    const double getRadiusOfGyration() const {
        return radiusOfGyration;
    }

    const Matrix3d & getTensor() const {
        return tensor;
    }

    const Vector3d & getTensorSingularValues() const {
        return tensorSingularValues;
    }

private:
    /**
     * Particle corodinates
     */
    const CoordMatrix & coords;

    /**
     * Gyration tensor
     */
    Matrix3d tensor;

    /**
     * Calculate gyration tensor and principle moments
     */
    void calculateGyrationTensorShapeDescriptors();

    /**
     * Singluar values/principal moments
     */
    Vector3d tensorSingularValues;

    /**
     * Shape descriptors
     */
    double radiusOfGyration, asphericity, acylindricity, anisotrophy;

};

} /* namespace Gape */

#endif /* SRC_UTIL_SHAPEMETRICS_H_ */
