//
// Created by gareth on 1/19/23.
//

#pragma once

#include <Geometry/point.h>
#include <string>

namespace Gape
{
	enum Geometry
	{
		Sphere,
		Vector,
		Arc,
		Cone
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
		VectorPharmFeatureGeometry(const RDGeom::Point3D& point1, const RDGeom::Point3D& point2);

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
	};

	/**
	 * Represents a feature that is a arc- typically an sp2 acceptor with two lone
	 * pairs.
	 * 
	 */
	class  ArcFeatureGeometry : public PharmFeatureGeometry
	{
		RDGeom::Point3D* points;

	public:

		/**
		 * The feature center is at point1 and point2 and point3 define the extent
		 * of the arc.
		 * 
		 * @param point1
		 * @param point2
		 * @param point3
		 */
		ArcFeatureGeometry(const RDGeom::Point3D& point1, const RDGeom::Point3D& point2, const RDGeom::Point3D& point3);

		std::string summary() const override;

		int getNumPoints() const override
		{
			return 3;
		}

		const RDGeom::Point3D& getPoint(int num) const override
		{
			return points[num];
		}

		~ArcFeatureGeometry()
		{
			delete[] points;
		}
	};

	/**
	 * Represents a feature that is a solid cone- typically an sp3 acceptor with
	 * three lone pairs.
	 * 
	 */
	class  ConeFeatureGeometry : public PharmFeatureGeometry
	{
		RDGeom::Point3D* points;

	public:

		/**
		 * point1 is the center of the cone and point2 and point3 are two points on
		 * the outside edge of the cone base such that the center of the circular
		 * base of the cone lies between points 2 and 3.
		 * 
		 * @param point1
		 * @param point2
		 * @param point3
		 */
		ConeFeatureGeometry(const RDGeom::Point3D& point1, const RDGeom::Point3D& point2, const RDGeom::Point3D& point3);

		std::string summary() const override;

		int getNumPoints() const override
		{
			return 3;
		}

		const RDGeom::Point3D& getPoint(int num) const override
		{
			return points[num];
		}

		~ConeFeatureGeometry()
		{
			delete[] points;
		}


	};


} // namespace Gape
