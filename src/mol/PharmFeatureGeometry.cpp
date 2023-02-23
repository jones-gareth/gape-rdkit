//
// Created by gareth on 1/19/23.
//

#include "PharmFeatureGeometry.h"
#include "../util/Util.h"

namespace Gape
{
	VectorPharmFeatureGeometry::VectorPharmFeatureGeometry(RDGeom::Point3D point1, RDGeom::Point3D point2):
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
} // Gape
