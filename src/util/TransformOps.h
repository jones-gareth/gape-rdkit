//
// Created by jones on 12/23/2022.
//

#ifndef GAPE_TRANSFORMOPS_H
#define GAPE_TRANSFORMOPS_H

#include <Geometry/point.h>
#include <Geometry/Transform3D.h>

namespace Gape
{
	using namespace RDKit;

	/**
	 * Determines the matrix transformation to perform a rotation (of angle
	 * angle) about an arbitrary axis (defined by p1 and p2). The result is
	 * returned in matrix.
	 *
	 * Ref: Newman & Sproull, Principles of Interactive computer Graphics,
	 * McGraw Hill, 1981 pp 346-348.
	 */
	void determineRotation(const RDGeom::Point3D& point1, const RDGeom::Point3D& point2, double angle,
	                       RDGeom::Transform3D& transform);

	/**
	 * Sets up two transformation matrices in and out. In moves point p1 to
	 * origin and aligns point p2 on z-axis. Out is the reverse
	 * transformation.
	 * 
	 * @param point1
	 * @param point2
	 * @param transformIn
	 * @param transformOut
	 */
	void setUpZAxisTransformation(const RDGeom::Point3D& point1, const RDGeom::Point3D& point2,
	                              RDGeom::Transform3D& transformIn, RDGeom::Transform3D& transformOut);


	/**
	 * Determines the angle between vectors (p1-p2) and (p3-p4).
	 * 
	 * @param p1
	 * @param p2
	 * @param p3
	 * @param p4
	 * @return
	 */
	double angleBetween(const RDGeom::Point3D& p1, const RDGeom::Point3D& p2, const RDGeom::Point3D& p3, const RDGeom::Point3D& p4);
} //namespace Gape
#endif //GAPE_TRANSFORMOPS_H
