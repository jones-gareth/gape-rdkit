/*
 * MolGeometry.cpp
 *
 *  Created on: Dec 16, 2015
 *      Author: Gareth Jones
 */

#include "MolGeometry.h"
#include "../util/CoordOps.h"

namespace GarethMol {

CoordVector generateOnePointForLinear(const Molecule & molecule,
		const Atom & atom1, const Atom & atom2, const double distance) {

	const auto coord1 = molecule.getCoord(atom1.getAtomNo());
	const auto coord2 = molecule.getCoord(atom2.getAtomNo());

	CoordVector diff = coord1 - coord2;
	diff = (distance * diff) / diff.norm();
	return coord1 + diff;
}

vector<CoordVector> generateThreePointsForTetrahedral(const Molecule & molecule,
		const Atom & atom1, const Atom & atom2, const double distance) {

	Transform<double, 3, Affine> in, out;
	transformToNegativeZaxis(molecule.getCoord(atom1.getAtomNo()),
			molecule.getCoord(atom2.getAtomNo()), in, out);

	// tetrahedral angle is 109.47
	auto angle = (19.47 * M_PI) / 180.0;
	auto z = distance * sin(angle);
	auto xy = distance * cos(angle);
	auto xyAngle = (30.0 * M_PI) / 180.0;
	auto x2 = xy * sin(xyAngle);
	auto y = -xy * cos(xyAngle);

	vector<CoordVector> rtn;
	rtn.reserve(3);
	rtn.push_back(static_cast<CoordVector>(in * CoordVector(xy, 0.0, z, 1.0)));
	rtn.push_back(static_cast<CoordVector>(in * CoordVector(x2, y, z, 1.0)));
	rtn.push_back(static_cast<CoordVector>(in * CoordVector(-x2, y, z, 1.0)));
	return rtn;
}

vector<CoordVector> generateTwoPointsForTetrahedral(const Molecule & molecule,
		const Atom & atom1, const Atom & atom2, const Atom & atom3,
		const double distance) {
	Transform<double, 3, Affine> in, out;
	// set up the transformation that moves atom1 to the origin and atom2 to the negative z axis
	transformToNegativeZaxis(molecule.getCoord(atom1.getAtomNo()),
			molecule.getCoord(atom2.getAtomNo()), in, out);
	// Find atom3 under this transformation
	CoordVector atom3Trans = in * molecule.getCoord(atom3.getAtomNo());
	atom3Trans.head(3).normalize();

	// spin atom3 round the z-azis by +/-120 degrees to create the two new points, then move back
	Transform<double, 3, Affine> zRot(
			AngleAxisd(2.0 * M_PI / 3.0, Vector3d::UnitZ()));
	Vector4d newPoint = zRot * atom3Trans;
	newPoint *= distance;
	newPoint[3] = 1.0;

	vector<CoordVector> rtn;
	rtn.reserve(2);
	rtn.push_back(static_cast<CoordVector>(out * newPoint));
	newPoint = zRot * newPoint;
	rtn.push_back(static_cast<CoordVector>(out * newPoint));
	return rtn;
}

CoordVector generateOnePointForTetrahedral(const Molecule & molecule,
		const Atom & atom1, const Atom & atom2, const Atom & atom3,
		const Atom & atom4, const double distance) {

	// find normalized vectors from vertices to center atom
	auto center = molecule.getCoord(atom1.getAtomNo()).head(3);
	auto point1 = molecule.getCoord(atom2.getAtomNo()).head(3);
	auto point2 = molecule.getCoord(atom3.getAtomNo()).head(3);
	auto point3 = molecule.getCoord(atom4.getAtomNo()).head(3);
	Vector3d vec1 = (center - point1);
	Vector3d vec2 = (center - point2);
	Vector3d vec3 = (center - point3);
	vec1.normalize();
	vec2.normalize();
	vec3.normalize();
	// determine centroid of these vectors
	Vector3d mean = (vec1 + vec2 + vec3) / 3.0;
	// find 4th vertex
	mean.normalize();
	mean = mean * distance;
	Vector3d newPoint = center + mean;
	return CoordVector(newPoint[0], newPoint[1], newPoint[2], 1.0);
}

CoordVector generateOnePointForTrigonal(const Molecule & molecule,
		const Atom & atom1, const Atom & atom2, const Atom & atom3,
		const double distance) {
	Transform<double, 3, Affine> in, out;
	// set up the transformation that moves atom1 to the origin and atom2 to the negative z axis
	transformToNegativeZaxis(molecule.getCoord(atom1.getAtomNo()),
			molecule.getCoord(atom2.getAtomNo()), in, out);
	// Find atom3 under this transformation
	Vector3d atom3Coord = molecule.getCoord(atom3.getAtomNo()).head(3);
	Vector3d atom3Trans = in * atom3Coord;
	atom3Trans.normalize();
	atom3Trans = atom3Trans * distance;
	// spin third atom round z-azis to generate new point
	Transform<double, 3, Affine> zRot(AngleAxisd(M_PI, Vector3d::UnitZ()));
	Vector3d newPoint = out * zRot * atom3Trans;
	return CoordVector(newPoint[0], newPoint[1], newPoint[2], 1.0);
}

vector<CoordVector> generateTwoPointsForTrigonal(const Molecule & molecule,
		const Atom & atom1, const Atom & atom2, const Atom & atom3,
		const double distance) {
	Transform<double, 3, Affine> in, out;
	// set up the transformation that moves atom1 to the origin and atom2 to the
	// negative z axis and atom3 to zx plane
	transformToZXplane(molecule.getCoord(atom1.getAtomNo()),
			molecule.getCoord(atom2.getAtomNo()), molecule.getCoord(atom3.getAtomNo()),in, out);
	// create new points
	auto thirtyDeg = M_PI / 6.0;
	auto x = distance * cos(thirtyDeg);
	auto z = distance * sin(thirtyDeg);
	CoordVector newPoint = CoordVector(x, 0.0, z, 1.0);
	vector<CoordVector> rtn;
	rtn.reserve(2);
	// move back to original corodinates
	newPoint = out * newPoint;
	rtn.push_back(newPoint);
	newPoint = CoordVector(-x, 0.0, z, 1.0);
	newPoint = out * newPoint;
	rtn.push_back(newPoint);
	return rtn;
}

vector<CoordVector> generateTwoPointsForIncompleteTrigonal(
		const Molecule & molecule, const Atom & atom1, const Atom & atom2,
		const double distance) {
	Transform<double, 3, Affine> in, out;
	// set up the transformation that moves atom1 to the origin and atom2 to the negative z axis
	transformToNegativeZaxis(molecule.getCoord(atom1.getAtomNo()),
			molecule.getCoord(atom2.getAtomNo()), in, out);
	// create new points on xz plane
	auto thirtyDeg = M_PI / 6.0;
	auto x = distance * cos(thirtyDeg);
	auto z = distance * sin(thirtyDeg);
	CoordVector newPoint = CoordVector(x, 0.0, z, 1.0);
	vector<CoordVector> rtn;
	rtn.reserve(2);
	// move back to original corodinates
	newPoint = out * newPoint;
	rtn.push_back(newPoint);
	newPoint = CoordVector(-x, 0.0, z, 1.0);
	newPoint = out * newPoint;
	rtn.push_back(newPoint);
	return rtn;
}

double outOfPlaneHeight(const Molecule & molecule, const Atom & atom,
		const Atom & neighbour1, const Atom & neighbour2,
		const Atom &neighbour3) {

    return GarethUtil::outOfPlaneHeight( molecule.getCoord(atom.getAtomNo()),
            molecule.getCoord(neighbour1.getAtomNo()),
            molecule.getCoord(neighbour2.getAtomNo()), molecule.getCoord(neighbour3.getAtomNo()));
}

}
