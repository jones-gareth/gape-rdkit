/*
 * CoordOps.cpp
 *
 *  Created on: May 20, 2014
 *      Author: Gareth Jones
 */

#include "CoordOps.h"
#include "Reporter.h"
#include "Util.h"
#include <cmath>

namespace GarethUtil {

using namespace Eigen;

bool checkSqrDistance(const double val, const CoordVector & a,
        const CoordVector & b) {
    double check = (a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1])
            + (a[2] - b[2]) * (a[2] - b[2]);
    return equals(val, check);
}

double sqrDistance(const CoordVector & a, const CoordVector & b) {

    Vector3d diff = a.head(3) - b.head(3);
    double sqrDistance = diff.squaredNorm();

    assert(checkSqrDistance(sqrDistance, a, b));
    return sqrDistance;
}

double angle(const CoordVector point1, CoordVector centerPoint,
        CoordVector point2) {
    Vector3d v1 = point1.head(3) - centerPoint.head(3);
    Vector3d v2 = point2.head(3) - centerPoint.head(3);
    auto dotProd = v1.dot(v2);
    double cosAngle = dotProd / (v1.norm() * v2.norm());
    return acos(cosAngle);
}

double distance(const CoordVector & a, const CoordVector & b) {
    double d = sqrt(sqrDistance(a, b));
    return d;
}

double rmsDistance(const CoordMatrix & a, const CoordMatrix & b) {
    auto diff = a - b;
    REPORT(Reporter::TRACE) << "Diff" << endl << diff;
    auto sqrDiff = diff.cwiseProduct(diff);
    auto sum = sqrDiff.sum();
    auto rmsd = sqrt(sum / static_cast<double>(diff.cols()));
    REPORT(Reporter::TRACE) << "sqrDiff " << endl << sqrDiff << endl
            << " rmsd " << rmsd;
    return rmsd;
}

const CoordVector getCoord(CoordMatrix & matrix, int index) {
    return matrix.col(index);
}

void transformToNegativeZaxis(const CoordVector & point1,
        const CoordVector & point2, Transform<double, 3, Affine> &in,
        Transform<double, 3, Affine> &out) {

    /* see Ref:
     Newman & Sproull, Principles of Interactive computer Graphics,
     McGraw Hill, 1981 pp 346-348. */

    // moves point1 to origin
    auto inTrans = Translation<double, 3>(-point1.head(3));
    // moves origin to point1
    auto outTrans = Translation<double, 3>(point1.head(3));

    // determine position of point 2 once point1 is at the origin
    Eigen::Vector3d diff = point2.head(3) - point1.head(3);
    diff.normalize();
    double a = diff(0);
    double b = diff(1);
    double c = diff(2);
    auto v = sqrt(b * b + c * c);

    double xAngle, yAngle;
    if (v == 0) {
        // point2 lies on x axis
        xAngle = .0;
        yAngle = a > 0 ? M_PI / 2.0 : -M_PI / 2.0;
    } else {
        // angle to rotate around x-axis to bring point 2 into zx plane
        xAngle = asin(b / v);
        if (c < 0)
            xAngle = -xAngle;
        // angle to rotate round y-axis to bring point 2 onto z axis
        yAngle = -atan(a / v);
        if (c >= 0)
            yAngle += M_PI;
        else
            yAngle = -yAngle;

#ifndef NDEBUG
        REPORT(Reporter::TRACE) << "a " << a << " b " << b << " c " << c;
        Transform<double, 3, Affine> xRot(
                AngleAxisd(xAngle, Vector3d::UnitX()));
        Transform<double, 3, Affine> yRot(
                AngleAxisd(yAngle, Vector3d::UnitY()));
        CoordVector newPoint = Transform<double, 3, Affine>(inTrans) * point2;
        REPORT(Reporter::TRACE) << "translated  " << point2 << " to "
                << newPoint;
        newPoint = xRot * newPoint;
        REPORT(Reporter::TRACE) << "rotated " << xAngle
                << " radians about x axis " << newPoint << " to " << newPoint;
        assert(equals(static_cast<double>(newPoint(1)), 0.0, 1e-10));
        newPoint = yRot * newPoint;
        REPORT(Reporter::TRACE) << "rotated " << yAngle
                << " radians about y axis to " << newPoint;
        assert(equals(static_cast<double>(newPoint(0)), 0.0, 1e-10));
        assert(equals(static_cast<double>(newPoint(1)), 0.0, 1e-10));
        assert(newPoint(2) < 0.0);
#endif

    }

    Matrix3d rot, invRot;
    rot = AngleAxisd(yAngle, Vector3d::UnitY())
            * AngleAxisd(xAngle, Vector3d::UnitX());
    invRot = AngleAxisd(-xAngle, Vector3d::UnitX())
            * AngleAxisd(-yAngle, Vector3d::UnitY());

    // build transforms into on x-axis and out
    in = rot * inTrans;
    out = outTrans * invRot;
}

void transformToZXplane(const CoordVector & point1, const CoordVector & point2,
        const CoordVector & point3, Transform<double, 3, Affine> &in,
        Transform<double, 3, Affine> &out) {
    // set up the transformation that moves point1 to the origin and point2 to the negative z axis
    transformToNegativeZaxis(point1, point2, in, out);
    // Find point under this transformation
    Vector3d point3Coord = point3.head(3);
    Vector3d point3Trans = in * point3Coord;
    // find rotation transform so that point3 now lies in zx plane (y-0)
    auto d = sqrt(
            static_cast<double>(point3Trans[0] * point3Trans[0]
                    + point3Trans[1] * point3Trans[1]));
    auto angle = acos(static_cast<double>(point3Trans[0]) / d);
    if (point3Trans[1] < 0) {
        angle = -angle;
    }

    // update out and in transforms with rotation to and from zx plane
    out = out * AngleAxisd(angle, Vector3d::UnitZ());
    in = AngleAxisd(-angle, Vector3d::UnitZ()) * in;

#ifndef NDEBUG
    point3Trans = in * point3Coord;
    double y = point3Trans[1];
    REPORT(Reporter::TRACE) << "y is " << y;
    assert(equals(y, 0, 1e-10));
    Vector3d test = out * point3Trans;
    for (int i = 0; i < 3; i++) {
        double c1 = point3Coord[i];
        double c2 = test[i];
        assert(equals(c1, c2, 1e-10));
    }
#endif

}

CoordMatrix centerCoordMatrix(const CoordMatrix &in) {
    double size = in.cols();
    CoordMatrix out;
    out.resize(4, in.cols());
    for (auto i = 0; i < 3; i++) {
        auto c = in.row(i).sum() / size;
        out.row(i) = in.row(i).array() - c;
    }
    out.row(3) = VectorXd::Constant(in.cols(), 1.0);
    return out;
}

CoordVector centroid(const CoordMatrix &in) {
    double size = in.cols();
    CoordVector out;
    for (auto i = 0; i < 3; i++) {
        auto c = in.row(i).sum() / size;
        out[i] = c;
    }
    out[3] = 1.0;
    return out;
}

double outOfPlaneHeight(const CoordVector & center,
        const CoordVector & neighbour1, const CoordVector & neighbour2,
        const CoordVector &neighbour3) {

    Transform<double, 3, Affine> in, out;
    // set up the transformation that moves atom1 to the origin and atom2 to the
    // negative z axis and atom3 to zx plane
    transformToZXplane(neighbour1, neighbour2, neighbour3, in, out);
    // transform center atom
    Vector3d atomCoord = center.head(3);
    Vector3d transAtom = in * atomCoord;
    // OOP distance is then just Y coordinate
    double y = transAtom[1];
    return fabs(y);
}

double torsionAngle(const CoordVector & a, const CoordVector & b,
        const CoordVector & c, const CoordVector & d) {

    Vector3d i = b.head(3) - a.head(3);
    Vector3d j = c.head(3) - b.head(3);
    Vector3d k = d.head(3) - c.head(3);

    // find triple vector product
    auto n1 = i.cross(j);
    auto n2 = j.cross(k);
    auto tvp = n1.cross(n2);

    auto sp = n2.dot(n1);
    auto sp2 = tvp.dot(j);
    auto mag = n1.norm() * n2.norm();
    auto cosVal = sp / mag;

    if (cosVal > 1) {
        cosVal = 1.0;
    } else if (cosVal < -1.0) {
        cosVal = -1.0;
    }

    auto angle = acos(cosVal);
    return sp2 < .0 ? -angle : angle;
}

double radiusOfGyration(const CoordMatrix & coords) {
    auto size = coords.cols();
    auto sum = .0;
    for (auto i = 0; i < size; i++) {
        for (auto j = 0; j < size; j++) {
            for (auto k =0; k<3; k++) {
                auto diff = coords(k, i) - coords(k, j);
                sum += diff*diff;
            }
        }
    }

    auto dSize = static_cast<double>(size);
    auto g = sum/(2.0*dSize*dSize);
    return g;
}

}
/* namespace GarethUtil */
