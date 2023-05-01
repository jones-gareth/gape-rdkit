//
// Created by gareth on 1/19/23.
//

#include "PharmFeatureGeometry.h"
#include "../util/Util.h"

namespace Gape
{
	VectorPharmFeatureGeometry::VectorPharmFeatureGeometry(const RDGeom::Point3D& point1,
	                                                       const RDGeom::Point3D& point2):
		PharmFeatureGeometry(Geometry::Vector)
	{
		points = new RDGeom::Point3D[2];
		points[0] = point1;
		points[1] = point2;
	}

	std::string VectorPharmFeatureGeometry::summary() const
	{
		return "{VECTOR " + toString(points[0]) + " " + toString(points[1])
			+ "}";
	}

	ArcFeatureGeometry::ArcFeatureGeometry(const RDGeom::Point3D& point1, const RDGeom::Point3D& point2,
	                                       const RDGeom::Point3D& point3) : PharmFeatureGeometry(Geometry::Arc)
	{
		points = new RDGeom::Point3D[3];
		points[0] = point1;
		points[1] = point2;
		points[2] = point3;
	}

	std::string ArcFeatureGeometry::summary() const
	{
		return "{ARC " + toString(points[0]) + " " + toString(points[1]) + " " + toString(points[2])
			+ "}";
	}

	ConeFeatureGeometry::ConeFeatureGeometry(const RDGeom::Point3D& point1, const RDGeom::Point3D& point2,
	                                         const RDGeom::Point3D& point3) : PharmFeatureGeometry(Geometry::Cone)
	{
		points = new RDGeom::Point3D[3];
		points[0] = point1;
		points[1] = point2;
		points[2] = point3;
	}

	std::string ConeFeatureGeometry::summary() const
	{
		return "{CONE " + toString(points[0]) + " " + toString(points[1]) + " " + toString(points[2])
			+ "}";
	}
} // Gape
