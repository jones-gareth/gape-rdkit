/*
 * CoordOps.h
 *
 *  Created on: May 20, 2014
 *      Author: Gareth Jones
 */

#ifndef COORDOPS_H_
#define COORDOPS_H_

#include <Eigen/Dense>

namespace GarethUtil {

using CoordMatrix = Eigen::Matrix<double, 4, Eigen::Dynamic, Eigen::ColMajor>;
using CoordVector = Eigen::Vector4d;
using namespace std;
using namespace Eigen;

//using CoordType = double;
//using CoordData = array<CoordType, 4>;

/**
 * Returns the squared distence between two vectors/coordinates
 * @param a
 * @param b
 * @return
 */
double sqrDistance(const CoordVector & a, const CoordVector & b);

/**
 * Return the distance between two coordinates
 *
 * @param a
 * @param b
 * @return
 */
double distance(const CoordVector & a, const CoordVector & b);

/**
 * Return the angle formed by three points.
 *
 * @param point1
 * @param centerPoint
 * @param point2
 * @return
 */
double angle(const CoordVector point1, CoordVector centerPoint, CoordVector point2);

/**
 * Return the RMS distance between two matrices
 *
 * @param a
 * @param b
 * @return
 */
double rmsDistance(const CoordMatrix & a, const CoordMatrix & b);


const CoordVector getCoord(CoordMatrix & matrix, int index);

/**
 * Creates transforms to move point1 onto the origin and point2 onto the negative z-axis
 *
 * @param point1
 * @param point2
 * @param in
 * @param out
 */
void transformToNegativeZaxis(const CoordVector & point1,
		const CoordVector & point2, Transform<double, 3, Affine> &in,
		Transform<double, 3, Affine> &out);

/**
 * Creates transforms to move point1 onto the origin, point2 onto the negative z-axis
 * and point 3 into the ZX plane (y = 0)
 *
 * @param point1
 * @param point2
 * @param point3
 * @param in
 * @param out
 */
void transformToZXplane(const CoordVector & point1, const CoordVector & point2,
		const CoordVector & point3, Transform<double, 3, Affine> &in,
		Transform<double, 3, Affine> &out);

/**
 * Centers the x,  y and z coordinates
 *
 * @param in
 * @return
 */
CoordMatrix centerCoordMatrix(const CoordMatrix &in);

/**
 * Returns the centroid of the coordinates
 *
 * @param in
 * @return
 */
CoordVector centroid(const CoordMatrix &in);

/**
 * Finds the distance center is from the plane formed by the three neighbours
 *
 * @param center
 * @param neighbour1
 * @param neighbour2
 * @param neighbour3
 * @return
 */
double outOfPlaneHeight( const CoordVector & center,
        const CoordVector & neighbour1, const CoordVector & neighbour2,
        const CoordVector &neighbour3);

/**
 * Returns the torsion angle betwwen 4 points.
 *
 * @param a
 * @param b
 * @param c
 * @param d
 * @return
 */
double torsionAngle(const CoordVector & a, const CoordVector & b,
        const CoordVector & c, const CoordVector & d);

/**
 * Returns the radius of gyration for the coordinates
 *
 * @param coords
 * @return
 */
double radiusOfGyration(const CoordMatrix & coords);


} /* namespace GarethUtil */

#endif /* COORDOPS_H_ */
