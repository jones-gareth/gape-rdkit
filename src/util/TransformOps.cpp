//
// Created by Gareth Jones on 12/23/2022.
//

#include "TransformOps.h"


#include <Geometry/point.h>

namespace Gape
{
	using namespace RDKit;

	void determineRotation(const RDGeom::Point3D& point1, const RDGeom::Point3D& point2, double angle,
	                       RDGeom::Transform3D& transform)
	{
		RDGeom::Point3D diff = point1 - point2;
		diff.normalize();
		RDGeom::Transform3D trans, xRot, yRot, rot, invTrans, invXRot, invYRot;
		trans.SetTranslation(-point1);
		invTrans.SetTranslation(point1);
		rot.SetRotation(angle, diff);
		transform.assign(invTrans * rot * trans);
	}


	/**
	 * Sets up two transformation matrices in and out. In moves point p1 to
	 * origin and alignes points p1 and p2 on z-axis. Out is the reverse
	 * transformation.
	 * 
	 * @param point1
	 * @param point2
	 * @param transformIn
	 * @param transformOut
	 */
	void setUpZAxisTransformation(const RDGeom::Point3D& point1, const RDGeom::Point3D& point2, RDGeom::Transform3D& transformIn, RDGeom::Transform3D& transformOut)
	{
		auto diff = point1 - point2;
		diff.normalize();
		const double a = diff.x;
		const double b = diff.y;
		const double c = diff.z;

		// trans: translates point one to the origin. xRot and yRot:
		// rotation about x and y axis such that point2 is aligned on
		// the z-axis. invTrans, invXrot and invYrot: inverse
		// translations and rotations.

		// Ref: Newman & Sproull, Principles of Interactive
		// computer Graphics, McGraw Hill, 1981 pp 346-348.

		RDGeom::Transform3D trans, invTrans, xRot, invXRot, yRot, invYRot;
		trans.SetTranslation(-point1);
		invTrans.SetTranslation(point1);
		const double v = sqrt(b * b + c * c);
		if (v > .0)
		{
			double* in = xRot.getData();
			double* out = invXRot.getData();
			in[5] = in[10] = out[5] = out[10] = c / v;
			in[6] = out[9] = b / v;
			in[9] = out[6] = -b / v;
		}
		double* in = yRot.getData();
		double* out = invYRot.getData();
		in[0] = in[10] = out[0] = out[10] = v;
		in[2] = out[8] = a;
		in[8] = out[2] = -a;

		transformIn.assign(yRot * xRot * trans);
		transformOut.assign(invTrans * invXRot * invYRot);
	}


	double angleBetween(const RDGeom::Point3D& p1, const RDGeom::Point3D& p2, const RDGeom::Point3D& p3, const RDGeom::Point3D& p4)
	{
		const auto v1 = p1 - p2;
		const auto v2 = p3 - p4;
		const auto angle = v1.angleTo(v2);
		return angle;
	}
}
