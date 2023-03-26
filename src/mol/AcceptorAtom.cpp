#include "AcceptorAtom.h"

#include "HydrogenBondingType.h"

namespace Gape
{
	AcceptorAtom::AcceptorAtom(const SuperpositionMolecule* spMol, const Atom* atom) : molecule(spMol), atom(atom)
	{
		auto& molAcceptors = molecule->getAcceptors();
		const auto it = molAcceptors.find(atom);
		assert(it != molAcceptors.end());
		auto hydrogenBondingType = it->second;

		int nLp = 0;
		const int degree = atom->getDegree();
		switch (atom->getHybridization())
		{
		case Atom::SP:
			nLp = 2 - degree;
			break;
		case Atom::SP2:
			nLp = 3 - degree;
			break;
		case Atom::SP3:
			nLp = 4 - degree;
			break;
		default:
			assert(false);
		}

		switch (hydrogenBondingType->geometry)
		{
		case HydrogenBondGeometry::None:
			nLp = 0;
			break;
		case HydrogenBondGeometry::Dir:
			break;
		
		}

		numberLonePairs = nLp;
	}

	void AcceptorAtom::addOnePairToLinear(const RDGeom::Point3D& atom1, const RDGeom::Point3D& atom2,
	                                      RDGeom::Point3D& lonePair, const double lonePairLength)
	{
		auto diff = atom1.directionVector(atom2);
		diff *= lonePairLength;
		lonePair = atom1 + diff;
	}

	void AcceptorAtom::addThreePairsToTetrahedral(const RDGeom::Point3D& origin, const RDGeom::Point3D& atom2,
	                                              RDGeom::Point3D& lonePair1, RDGeom::Point3D& lonePair2,
	                                              RDGeom::Point3D& lonePair3, double lonePairLength)
	{

	}
}
