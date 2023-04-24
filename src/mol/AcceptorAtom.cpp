#include "AcceptorAtom.h"

#include "HydrogenBondingType.h"
#include "../util/TransformOps.h"
#include "../util/Reporter.h"

#include <boost/format.hpp>

namespace Gape
{
	const double AcceptorAtom::lonePairLength = 2.9;

	AcceptorAtom::AcceptorAtom(const SuperpositionMolecule* spMol, const Atom* atom) : molecule(spMol), atom(atom)
	{
		auto& molAcceptors = molecule->getAcceptors();
		const auto it = molAcceptors.find(atom);
		assert(it != molAcceptors.end());
		hydrogenBondingType = it->second;

		numberLonePairs = countLonePairs();
		addLonePairs();
	}

	void AcceptorAtom::addOnePairToLinear(const RDGeom::Point3D& atom1, const RDGeom::Point3D& atom2,
	                                      RDGeom::Point3D& lonePair, const double lonePairLength)
	{
		auto diff = atom1.directionVector(atom2);
		diff *= lonePairLength;
		lonePair = atom1 + diff;
	}

	void AcceptorAtom::addThreePairsToTetrahedral(const RDGeom::Point3D& origin, const RDGeom::Point3D& atom1,
	                                              RDGeom::Point3D& lonePair1, RDGeom::Point3D& lonePair2,
	                                              RDGeom::Point3D& lonePair3, const double lonePairLength)
	{
		// Set up transformation that tranforms donor to origin and atom2
		// along z-axis
		RDGeom::Transform3D transformIn, transformOut, rotation;
		setUpZAxisTransformation(origin, atom1, transformIn, transformOut);

		const double xv = lonePairLength * cos((19.47 / 180 * M_PI));
		const double zv = lonePairLength * sin((19.47 / 180 * M_PI));
		RDGeom::Point3D dummy(xv, .0, zv);
		// generate 1st lone pair at 0 degrees
		rotation.SetRotation(0, RDGeom::Z_Axis);
		rotation.TransformPoint(dummy);
		lonePair1 = dummy;
		// rotate by 120 degrees to generate 2nd and 3rd lone pairs.
		rotation.SetRotation(2.0 * M_PI / 3.0, RDGeom::Z_Axis);
		rotation.TransformPoint(dummy);
		lonePair2 = dummy;
		rotation.TransformPoint(dummy);
		lonePair3 = dummy;

		// Apply inverse transformations to generate the final
		// lone pair co-ordinates
		transformOut.TransformPoint(lonePair1);
		transformOut.TransformPoint(lonePair2);
		transformOut.TransformPoint(lonePair1);
	}

	void AcceptorAtom::addTwoPairsToTetrahedral(const RDGeom::Point3D& origin, const RDGeom::Point3D& atom1,
	                                            const RDGeom::Point3D& atom2, RDGeom::Point3D& lonePair1,
	                                            RDGeom::Point3D& lonePair2, const double lonePairLength)
	{
		// set up the the transformation that moves the acceptor to the
		// origin and atom1 on the z-axis
		RDGeom::Transform3D transformIn, transformOut, rotation;
		setUpZAxisTransformation(origin, atom1, transformIn, transformOut);

		// dummy is the co-ordinates of atom2 under this transformation
		RDGeom::Point3D dummy(atom2);
		transformIn.TransformPoint(dummy);
		/* find the angle dummy makes with the zy-axis */
		double angle = asin(dummy[0] / (dummy[0] * dummy[0] + dummy[1] * dummy[1]));
		if (dummy[1] < 0.0)
			angle = M_PI - angle;
		angle *= -1.0;
		// dummy1 is dummy rotated onto the zy plane
		rotation.SetRotation(angle, RDGeom::Z_Axis);
		RDGeom::Point3D dummy1(dummy);
		rotation.TransformPoint(dummy1);
		// atom1 lies along -ve z-axis; atom2 should be in +ve quadrant
		// of zy plane

		// Alpha is the angle between atom1-acceptor-atom2
		const double alpha = M_PI / 2 + asin(dummy1[2] / sqrt(dummy1[1] * dummy1[1] + dummy1[2] * dummy1[2]));
		// Now add the two lone pairs in a plane perpendicular to the
		// zy plane that passes through the origin and bisects
		// alpha. Beta is the angle out between this plane and the y
		// axis.
		const double beta = 0.5 * (2.0 * M_PI - alpha) - 0.5 * M_PI;
		// The angle between the two lone pairs is 109.47 so we can determine
		// their projections in the x and zy planes, as one lone pair will be
		// added above the zy plane and one lone pair will be added an equal
		// distance below it
		const double xLen = lonePairLength * sin(54.735 * M_PI / 180);
		const double zyLen = lonePairLength * cos(54.735 * M_PI / 180);
		lonePair1[0] = xLen;
		lonePair2[0] = -xLen;
		// The projections of the lone pairs on the z and y axis are
		// determined from beta
		lonePair1[1] = lonePair2[1] = -1.0 * zyLen * cos(beta);
		lonePair1[2] = lonePair2[2] = 1.0 * zyLen * sin(beta);
		lonePair1[3] = lonePair2[3] = 1.0;

		/* Move lone pairs back to the original co-ordinate system. */
		rotation.SetRotation(-1.0 * angle, RDGeom::Z_Axis);
		rotation.TransformPoint(lonePair1);
		rotation.TransformPoint(lonePair2);

		transformOut.TransformPoint(lonePair1);
		transformOut.TransformPoint(lonePair2);
	}

	/**
	 * rescaleAtom2 is on the line between atom1 and atom2, unit distance from
	 * atom1.
	 *
	 * @param atom1
	 * @param atom2
	 * @param rescaleAtom2
	 */
	void resaleBond(const RDGeom::Point3D& atom1, const RDGeom::Point3D& atom2, RDGeom::Point3D& rescaleAtom2)
	{
		auto diff = atom2 - atom1;
		diff.normalize();
		rescaleAtom2 = atom1 + diff;
	}

	void AcceptorAtom::addOnePairToTetrahedral(const RDGeom::Point3D& origin, const RDGeom::Point3D& atom1,
	                                           const RDGeom::Point3D& atom2, const RDGeom::Point3D& atom3,
	                                           RDGeom::Point3D& lonePair1, const double lonePairLength)
	{
		// rescale bondLengths to 1
		RDGeom::Point3D rescaleAtom1, rescaleAtom2, rescaleAtom3;
		resaleBond(origin, atom1, rescaleAtom1);
		resaleBond(origin, atom2, rescaleAtom2);
		resaleBond(origin, atom3, rescaleAtom3);

		const double xv = (rescaleAtom1[0] + rescaleAtom2[0] + rescaleAtom3[0]) / 3.0;
		const double yv = (rescaleAtom1[1] + rescaleAtom2[1] + rescaleAtom3[1]) / 3.0;
		const double zv = (rescaleAtom1[2] + rescaleAtom2[2] + rescaleAtom3[2]) / 3.0;
		const RDGeom::Point3D centroid(xv, yv, zv);
		// find lone pair
		auto lpVec = origin - centroid;
		lpVec.normalize();
		lpVec = lpVec * lonePairLength;
		lonePair1 = origin + lpVec;
	}


	bool AcceptorAtom::addOnePairToTrigonal(const RDGeom::Point3D& origin, const RDGeom::Point3D& atom1,
	                                        const RDGeom::Point3D& atom2, RDGeom::Point3D& lonePair1,
	                                        const double lonePairLength)
	{
		// Set up the transformation that puts <origin> at the origin
		// and atom1 along the z-axis.
		RDGeom::Transform3D transformIn, transformOut, rotation;
		setUpZAxisTransformation(origin, atom1, transformIn, transformOut);

		// move atom2
		RDGeom::Point3D atom2Moved(atom2);
		transformIn.TransformPoint(atom2Moved);
		// <angle> is the angle atom2 makes with the zy plane */
		double angle = asin(atom2Moved[0]
			/ (sqrt(atom2Moved[0] * atom2Moved[0] + atom2Moved[1] * atom2Moved[1])));
		if (atom2Moved[1] < 0.0)
			angle = M_PI - angle;
		angle *= -1.0;
		// new atom2 is rotated to dummy1 which lies on the zy plane
		rotation.SetRotation(angle, RDGeom::Z_Axis);
		rotation.TransformPoint(atom2Moved);
		// atom1 lies along -ve z-axis; atom2 should be in +ve
		// quadrant of zy plane
		if (atom2Moved[1] < 0 || atom2Moved[2] < 0)
			return false;

		// alpha is the angle between atom1-origin-atom2
		const double alpha = M_PI / 2 + asin(
			atom2Moved[2] / sqrt(atom2Moved[1] * atom2Moved[1] + atom2Moved[2] * atom2Moved[2]));
		// beta is the angle we want the lone pair to make with the
		// z-axis
		const double beta = 0.5 * (2.0 * M_PI - alpha) - 0.5 * M_PI;
		lonePair1[0] = 0;
		lonePair1[1] = -1.0 * lonePairLength * cos(beta);
		lonePair1[2] = lonePairLength * sin(beta);

		// move lone pair back
		rotation.SetRotation(angle, RDGeom::Z_Axis);
		rotation.TransformPoint(lonePair1);
		transformOut.TransformPoint(lonePair1);
		return true;
	}

	bool AcceptorAtom::addTwoPairsToTrigonal(const RDGeom::Point3D& atom1, const RDGeom::Point3D& origin,
	                                         const RDGeom::Point3D& atom2, RDGeom::Point3D& lonePair1,
	                                         RDGeom::Point3D& lonePair2, const double lonePairLength)
	{
		// Set up the transformation that moves origin to the origin
		// and puts atom1 on the z-axis
		RDGeom::Transform3D transformIn, transformOut, rotation;
		setUpZAxisTransformation(origin, atom1, transformIn, transformOut);

		RDGeom::Point3D dummy(atom2), dummy2(atom1);
		// move atom1 and atom2 to dummy2 and dummy */
		transformIn.TransformPoint(dummy);
		transformIn.TransformPoint(dummy2);

		// angle is the angle atom2 makes with the zy plane
		double angle = asin(dummy[0] / (sqrt(dummy[0] * dummy[0] + dummy[1] * dummy[1])));
		if (dummy[1] < 0.0)
			angle = M_PI - angle;
		angle *= -1.0;
		rotation.SetRotation(angle, RDGeom::Z_Axis);
		rotation.TransformPoint(dummy);

		// atom1 lies along -ve z-axis; atom2 should be in +ve
		// quadrant of zy plane
		if (dummy[1] < 0 || dummy[2] < 0)
			return false;

		// can now add lone pairs onto atom1 in the zy plane so that
		// the angle origin-atom1-lp is 120 degrees
		constexpr double alpha = M_PI / 3.0;
		const double l = lonePairLength * sin(alpha);
		lonePair1[0] = lonePair2[0] = 0;
		lonePair1[1] = l;
		lonePair2[1] = -l;
		lonePair1[2] = lonePair2[2] = dummy2[2] - lonePairLength * cos(alpha);
		lonePair1[3] = lonePair2[3] = 1.0;

		// move lone pairs back
		rotation.SetRotation(-1.0 * angle, RDGeom::Z_Axis);
		rotation.TransformPoint(lonePair1);
		rotation.TransformPoint(lonePair2);
		transformOut.TransformPoint(lonePair1);
		transformOut.TransformPoint(lonePair2);
		return true;
	}

	void AcceptorAtom::addTwoPairsRandomlyToTrigonal(const RDGeom::Point3D& atom1, const RDGeom::Point3D& origin,
	                                                 RDGeom::Point3D& lonePair1, RDGeom::Point3D& lonePair2,
	                                                 const double alpha, const double lonePairLength)
	{
		// get the transformation that moves <origin> to the origin
		// and atom1 to the z-axis
		RDGeom::Transform3D transformIn, transformOut, rotation;
		setUpZAxisTransformation(origin, atom1, transformIn, transformOut);

		// Dummy is atom1 transformed
		RDGeom::Point3D dummy(atom1);
		transformIn.TransformPoint(dummy);
		// Atom1 lies along -ve z-axis. Add the lone pairs in the zy
		// plane to atom1 such that the angle lone-pair-atom1-origin
		// is 120 degrees
		const double l = lonePairLength * sin(alpha);
		lonePair1[0] = lonePair2[0] = 0;
		lonePair1[1] = l;
		lonePair2[1] = -l;
		lonePair1[2] = lonePair2[2] = dummy[2] - lonePairLength
			* cos(alpha);
		lonePair1[3] = lonePair2[3] = 1.0;

		// Angle could be random- but we just make it 0
		double angle = 0;
		rotation.SetRotation(angle, RDGeom::Z_Axis);
		rotation.TransformPoint(lonePair1);
		rotation.TransformPoint(lonePair2);
		// move lone pairs back
		transformOut.TransformPoint(lonePair1);
		transformOut.TransformPoint(lonePair2);
	}

	int AcceptorAtom::countTrueLonePairs() const
	{
		const int nConnections = atom->getDegree();
		int total = 0;

		// hardwire Nitro oxygens
		if (molecule->isNitroOxygen(*atom))
		{
			total = 3;
		}
		else if (atom->getHybridization() == Atom::SP)
		{
			total = 2;
		}
		else if (atom->getHybridization() == Atom::SP2)
		{
			total = 3;
		}
		else if (atom->getHybridization() == Atom::SP3)
		{
			total = 4;
		}

		return total - nConnections;
	}

	int AcceptorAtom::countLonePairs() const
	{
		switch (hydrogenBondingType->geometry)
		{
		case Dir:
			return countTrueLonePairs();
		case Plane:
			{
				const auto numLP = countTrueLonePairs();
				if (numLP != 2)
				{
					throw std::domain_error(
						(boost::format("Can't create plane lone pair geometry on acceptor atom %d") % atom->getIdx()).
						str());
				}
				return 2;
			}
		case Cone:
		case None:
			// For the cone and none the lone pair is just a
			// representation of forward geometry and not a real
			// lone pair.
			return 1;
		}

		// shouldn't get here
		return 0;
	}

	void AcceptorAtom::add1LpToSp3()
	{
		if (auto degree = atom->getDegree(); degree != 3)
		{
			throw std::domain_error(
				(boost::format("add1LpToSp3: atom %d has %d connections") % atom->getIdx() % degree).str());
		}

		const auto& conformer = molecule->getReferenceConformer();
		const auto& coord = conformer.getAtomPos(atom->getIdx());
		std::vector<RDGeom::Point3D> otherCoords;
		for (const auto nbr : molecule->getMol().atomNeighbors(atom))
		{
			otherCoords.push_back(conformer.getAtomPos(nbr->getIdx()));
		}

		RDGeom::Point3D lonePair;
		addOnePairToTetrahedral(coord, otherCoords[0], otherCoords[1], otherCoords[2], lonePair, lonePairLength);
		lonePairs.clear();
		lonePairs.push_back(lonePair);
	}

	void AcceptorAtom::add1LpToSp2()
	{
		if (auto degree = atom->getDegree(); degree != 2)
		{
			throw std::domain_error(
				(boost::format("add1LpToSp2: atom %d has %d connections") % atom->getIdx() % degree).str());
		}

		const auto& conformer = molecule->getReferenceConformer();
		const auto& coord = conformer.getAtomPos(atom->getIdx());
		std::vector<RDGeom::Point3D> otherCoords;
		for (const auto nbr : molecule->getMol().atomNeighbors(atom))
		{
			otherCoords.push_back(conformer.getAtomPos(nbr->getIdx()));
		}

		RDGeom::Point3D lonePair;
		addOnePairToTrigonal(coord, otherCoords[0], otherCoords[1], lonePair, lonePairLength);
		lonePairs.clear();
		lonePairs.push_back(lonePair);
	}

	void AcceptorAtom::addLpToSp1()
	{
		if (auto degree = atom->getDegree(); degree != 1)
		{
			throw std::domain_error(
				(boost::format("add1LpToSp1: atom %d has %d connections") % atom->getIdx() % degree).str());
		}

		const auto& conformer = molecule->getReferenceConformer();
		const auto& coord = conformer.getAtomPos(atom->getIdx());
		const auto& atom2 = molecule->getMol().atomNeighbors(atom).begin().current;
		const auto& coord2 = conformer.getAtomPos(atom2->getIdx());

		RDGeom::Point3D lonePair;
		addOnePairToLinear(coord, coord2, lonePair, lonePairLength);
		lonePairs.clear();
		lonePairs.push_back(lonePair);
	}

	void AcceptorAtom::add2LpToSp3()
	{
		if (auto degree = atom->getDegree(); degree != 2)
		{
			throw std::domain_error(
				(boost::format("add2LpToSp3: atom %d has %d connections") % atom->getIdx() % degree).str());
		}

		const auto& conformer = molecule->getReferenceConformer();
		const auto& coord = conformer.getAtomPos(atom->getIdx());
		std::vector<RDGeom::Point3D> otherCoords;
		for (const auto nbr : molecule->getMol().atomNeighbors(atom))
		{
			otherCoords.push_back(conformer.getAtomPos(nbr->getIdx()));
		}

		RDGeom::Point3D lonePair1, lonePair2;
		addTwoPairsToTetrahedral(coord, otherCoords[0], otherCoords[1], lonePair1, lonePair2, lonePairLength);
		lonePairs.clear();
		lonePairs.push_back(lonePair1);
		lonePairs.push_back(lonePair2);
	}

	void AcceptorAtom::add2LpToSp2()
	{
		if (auto degree = atom->getDegree(); degree != 1)
		{
			throw std::domain_error(
				(boost::format("add2LpToSp2: atom %d has %d connections") % atom->getIdx() % degree).str());
		}

		const auto& conformer = molecule->getReferenceConformer();
		const auto& mol = molecule->getMol();
		const auto otherAtom = mol.atomNeighbors(atom).begin().current;
		const Atom* thirdAtom = nullptr;
		for (const auto a : mol.atomNeighbors(otherAtom))
		{
			if (a->getIdx() != atom->getIdx())
			{
				thirdAtom = a;
			}
		}
		assert(thirdAtom != nullptr);

		const auto& coord = conformer.getAtomPos(atom->getIdx());
		const auto& coord2 = conformer.getAtomPos(otherAtom->getIdx());
		const auto& coord3 = conformer.getAtomPos(thirdAtom->getIdx());
		RDGeom::Point3D lonePair1, lonePair2;
		if (!addTwoPairsToTrigonal(coord, coord2, coord3, lonePair1, lonePair2, lonePairLength))
		{
			throw std::domain_error(
				(boost::format("add2LpToSp2: failed to add two pairs to trigonal for atom %d") % atom->getIdx()).str()
			);
		}
		lonePairs.clear();
		lonePairs.push_back(lonePair1);
		lonePairs.push_back(lonePair2);
	}

	void AcceptorAtom::addLpToOco2()
	{
		if (auto degree = atom->getDegree(); degree != 1)
		{
			throw std::domain_error(
				(boost::format("add1LpToOco2: atom %d has %d connections") % atom->getIdx() % degree).str());
		}

		const auto& conformer = molecule->getReferenceConformer();
		const auto& mol = molecule->getMol();
		const auto otherAtom = mol.atomNeighbors(atom).begin().current;
		const Atom* thirdAtom = nullptr;
		for (const auto a : mol.atomNeighbors(otherAtom))
		{
			if (a->getIdx() != atom->getIdx() && a->getAtomicNum() == 8 && a->getDegree() == 1)
			{
				thirdAtom = a;
			}
		}
		if (thirdAtom == nullptr)
		{
			throw std::domain_error(
				(boost::format("addLpToOco2: can't find other oco2 paired to atom %d") % atom->getIdx()).str()
			);
		}

		const auto& coord = conformer.getAtomPos(atom->getIdx());
		const auto& coord2 = conformer.getAtomPos(otherAtom->getIdx());
		const auto& coord3 = conformer.getAtomPos(thirdAtom->getIdx());
		RDGeom::Point3D lonePair1, lonePair2;
		if (!addTwoPairsToTrigonal(coord, coord2, coord3, lonePair1, lonePair2, lonePairLength))
		{
			throw std::domain_error(
				(boost::format("add2LpToOco2: failed to add two pairs to trigonal for atom %d") % atom->getIdx()).str()
			);
		}
		lonePairs.clear();
		lonePairs.push_back(lonePair1);
		lonePairs.push_back(lonePair2);
	}

	void AcceptorAtom::add3LpToSp3()
	{
		if (auto degree = atom->getDegree(); degree != 1)
		{
			throw std::domain_error(
				(boost::format("add3LpToSp3: atom % d has % d connections") % atom->getIdx() % degree).str());
		}

		const auto& conformer = molecule->getReferenceConformer();
		const auto& coord = conformer.getAtomPos(atom->getIdx());
		const auto& atom2 = molecule->getMol().atomNeighbors(atom).begin().current;
		const auto& coord2 = conformer.getAtomPos(atom2->getIdx());

		RDGeom::Point3D lonePair1, lonePair2, lonePair3;
		addThreePairsToTetrahedral(coord, coord2, lonePair1, lonePair2, lonePair3, lonePairLength);
		lonePairs.clear();
		lonePairs.push_back(lonePair1);
		lonePairs.push_back(lonePair2);
		lonePairs.push_back(lonePair3);
	}

	void AcceptorAtom::addAtomLonePairs()
	{
		const auto geometry = hydrogenBondingType->geometry;

		// Hardwire carboxylic acid
		if (molecule->isCarboxylateOxygen(*atom))
		{
			addLpToOco2();
		}
		// hardwire nitro group
		else if (molecule->isNitroOxygen(*atom))
		{
			addLpToOco2();
		}
		else
		{
			if (atom->getHybridization() == Atom::SP3)
			{
				if (numberLonePairs == 1)
					add1LpToSp3();
				else if (numberLonePairs == 2)
					add2LpToSp3();
				else if (numberLonePairs == 3)
					add3LpToSp3();
			}
			else if (atom->getHybridization() == Atom::SP2)
			{
				if (numberLonePairs == 1)
					add1LpToSp2();
				else if (numberLonePairs == 2)
					add2LpToSp2();
			}
			else if (atom->getHybridization() == Atom::SP)
			{
				if (numberLonePairs == 1)
					addLpToSp1();
			}
			else
			{
				throw std::domain_error(
					(boost::format("addAtomLonePairs: Unable to add lone pairs to atom %d") % atom->getIdx()).str()
				);
			}
		}

		assert(checkLonePairLengths());
	}

	void AcceptorAtom::addOnePairToAtom()
	{
		switch (atom->getDegree())
		{
		case 1:
			addLpToSp1();
			break;
		case 2:
			add1LpToSp2();
			break;
		case 3:
			add1LpToSp3();
			break;
		default:
			throw std::domain_error(
				(boost::format("addOnePairToAtom: Unable to add a single LP to atom %d with degree %d") % atom->getIdx()
					% atom->getDegree()).str()
			);
		}
	}

	void AcceptorAtom::addLonePairs()
	{
		switch (hydrogenBondingType->geometry)
		{
		case Dir:
			addAtomLonePairs();
			break;
		case Plane:
			assert(numberLonePairs == 2);
			addAtomLonePairs();
			break;
		case Cone:
		case None:
			assert(numberLonePairs == 1);
			addOnePairToAtom();
			break;
		}
	}

	bool AcceptorAtom::checkLonePairLengths()
	{
		const auto& atomCoord = molecule->getReferenceConformer().getAtomPos(atom->getIdx());
		for (const auto& lonePair: lonePairs)
		{
			const auto diff = lonePair - atomCoord;
			const auto len = diff.length();
			if (abs(lonePairLength - len) > 0.01)
			{
				REPORT(Reporter::WARN) << "Lone pair length error: expected " << lonePairLength << " got " << len;
				return false;
			}
		}

		if (numberLonePairs != lonePairs.size())
		{
			REPORT(Reporter::WARN) << "Lone pair count error: expected " << numberLonePairs << " got " << lonePairs.size();
			return false;
		}

		return true;
	}

} // namespace Gape
