//
// Created by gareth on 1/19/23.
//

#ifndef GAPE_PHARMFEATUREGEOMETRY_H
#define GAPE_PHARMFEATUREGEOMETRY_H

#include <Geometry/point.h>
#include <string>

namespace Gape
{
	enum Geometry
	{
		Sphere,
		Vector,
		Arc,
		MultiVector
	};

	class PharmFeatureGeometry
	{
	protected:
		const Geometry geometry;

	public:
		PharmFeatureGeometry(Geometry g) : geometry(g)
		{
		}

		virtual std::string summary() const = 0;
		/**
		 * Gets the number of points used by this feature.
		 *
		 * @return
		 */
		virtual int getNumPoints() const = 0;

		/**
		 * Returns the coordinates of a point
		 *
		 * @param num
		 * @return
		 */
		virtual const RDGeom::Point3D& getPoint(int num) const = 0;
	};


	/**
	 * Represents a pharmacophore feature that is a vector. For example, a donor
	 * hydrogen or an sp1 acceptor.
	 *
	 */
	class VectorPharmFeatureGeometry : public PharmFeatureGeometry
	{
		RDGeom::Point3D* points;

	public:
		/**
		 * Constructor. The feature directionality is from point 1 to point 2.
		 *
		 * @param point1
		 * @param point2
		 */
		VectorPharmFeatureGeometry(RDGeom::Point3D point1, RDGeom::Point3D point2);

		std::string summary() const override;

		int getNumPoints() const override
		{
			return 2;
		}

		const RDGeom::Point3D& getPoint(int num) const override
		{
			return points[num];
		}

		~VectorPharmFeatureGeometry()
		{
			delete[] points;
		}
	}; // Gape
}
#endif //GAPE_PHARMFEATUREGEOMETRY_H
