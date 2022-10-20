/*
 * Atom.cpp
 *
 *  Created on: Apr 28, 2014
 *      Author: gjones
 */

#include <boost/format.hpp>

#include "Atom.h"
#include "Bond.h"
#include "Molecule.h"
#include "Ring.h"
#include "../util/Reporter.h"
#include "SubstructureSearch.h"
#include "SmartsParser.h"
#include "MolGeometry.h"

namespace GarethMol {

using namespace std;
using namespace GarethUtil;

using AtomTypeId = AtomType::AtomTypeId;
using BondTypeId = BondType::BondTypeId;

const string Atom::info() const {
	return to_string(atomNo + 1) + ":" + atomType->getName() + " [" + name + "]";
}

const bool Atom::sameAtomType(const Atom & otherAtom) const {
	return getAtomType().getType() == otherAtom.getAtomType().getType();
}

void AtomNeighbourhood::create(const Atom & atom, const Molecule & molecule) {
	auto atomNo = atom.getAtomNo();
	auto filter =
			[atomNo] (const Bond & bond) {
				return bond.getAtom1().getAtomNo() == atomNo || bond.getAtom2().getAtomNo() == atomNo;
			};
	bonds = filterUniquePtrListToConstRawPtrList<Bond>(molecule.getBonds(),
			filter);
	auto realBondFilter = [] (const Bond * bond) {
		return AtomType::isRealAtom(bond->getAtom1().getAtomTypeId()) &&
		AtomType::isRealAtom(bond->getAtom2().getAtomTypeId());
	};
	atomBonds = filterListToNewList<const Bond *>(bonds, realBondFilter);

	auto atomMapper = [atomNo] (const Bond * bond) -> const Atom * {
		if (bond->getAtom1().getAtomNo() == atomNo)
		return bond->getAtom2Ptr();
		return bond->getAtom1Ptr();
	};
	allNeighbours = mapToNewList<const Bond *, const Atom *>(bonds, atomMapper);

	auto heavyFilter = [] (const Atom * atom) {
		return AtomType::isHeavy(atom->getAtomTypeId());
	};
	heavyAtomNeighbours = filterListToNewList<const Atom *>(allNeighbours,
			heavyFilter);

	auto realFilter = [] (const Atom * atom) {
		return AtomType::isRealAtom(atom->getAtomTypeId());
	};
	atomNeighbours = filterListToNewList<const Atom *>(allNeighbours,
			realFilter);

	auto accum = [] (int count, const Atom * atom) {
		auto inc = atom->getAtomTypeId() == AtomTypeId::H ? 1 : 0;
		return count + inc;
	};
	noHydrogens = reduce<const Atom *, int>(atomNeighbours, 0, accum);

	auto angleFilter = [atomNo] (const Angle & angle) {
		return angle.getAtom1() == atomNo || angle.getAtom2() == atomNo
		|| angle.getAtom3() == atomNo;
	};
	angles = filterUniquePtrListToConstRawPtrList<Angle>(molecule.getAngles(),
				angleFilter);

	auto torsionFilter = [atomNo] (const Torsion & torsion) {
		return torsion.getAtom1() == atomNo || torsion.getAtom2() == atomNo
				|| torsion.getAtom3() == atomNo || torsion.getAtom4() == atomNo;
	};
	torsions = filterUniquePtrListToConstRawPtrList<Torsion>(molecule.getTorsions(),
				torsionFilter);
}

const int AtomNeighbourhood::countMetalNeighbours() const {
	auto accum = [] (int count, const Atom * atom) {
		auto inc = atom->getAtomType().isMetal() ? 1 : 0;
		return count + inc;
	};
	return reduce<const Atom *, int>(atomNeighbours, 0, accum);
}

const int bondCount(const vector<const Bond *> & bonds,
		const BondType::BondTypeId type) {
	auto accum = [type] (int count, const Bond * bond) -> int {
		auto inc = bond->getBondTypeId() == type ? 1 : 0;
		return count + inc;
	};
	return reduce<const Bond *, int>(bonds, 0, accum);
}

const int AtomNeighbourhood::countAromaticBonds() const {
	return bondCount(atomBonds, BondTypeId::AR);
}

const int AtomNeighbourhood::countSingleAndAmideBonds() const {
	auto accum =
			[] (int count, const Bond * bond) -> int {
				auto bondType = bond->getBondTypeId();
				auto inc = (bondType == BondTypeId::SINGLE || bondType == BondTypeId::AM) ? 1 : 0;
				return count + inc;
			};
	return reduce<const Bond *, int>(atomBonds, 0, accum);
}

const int AtomNeighbourhood::countDoubleBonds() const {
	return bondCount(atomBonds, BondTypeId::DOUBLE);
}

const int AtomNeighbourhood::countTripleBonds() const {
	return bondCount(atomBonds, BondTypeId::TRIPLE);
}

const int AtomNeighbourhood::totalBondOrder() const {
	auto accum = [] (double count, const Bond * bond) {
		return count + BondType::bondOrder(bond->getBondTypeId());
	};
	double bondOrder = reduce<const Bond *, double>(atomBonds, .0, accum);
	return round(bondOrder);
}

const int AtomNeighbourhood::countDegree(bool countImplicit) const {
	if (countImplicit) {
		return atomNeighbours.size();
	}
	auto accum = [] (int count, const Atom * atom) {
		auto inc = atom->isImplicit() ? 0 : 1;
		return count + inc;
	};
	return reduce<const Atom *, int>(atomNeighbours, 0, accum);
}

const int AtomNeighbourhood::countImplicitHydrogens() const {
	auto accum = [] (int count, const Atom * atom) {
		auto inc = atom->isImplicit() &&
		atom->getAtomTypeId() == AtomTypeId::H ? 1 : 0;
		return count + inc;
	};
	return reduce<const Atom *, int>(atomNeighbours, 0, accum);
}

const int Atom::totalBondOrder(const Molecule & molecule) const {
	return totalBondOrder(molecule.getBonds());
}

const int Atom::totalBondOrder(
		const std::vector<unique_ptr<Bond>> & bonds) const {
	struct OrderInfo {
		double order;
		int nAromatic;
	};
	OrderInfo orderInfo = { 0.0, 0 };

	auto accum =
			[this] (OrderInfo & info, const unique_ptr<Bond> & bond) {
				if (bond->getAtom1().getAtomNo()==atomNo|| bond->getAtom2().getAtomNo() == atomNo) {
					REPORT(Reporter::TRACE) << bond->info() << " order is " << BondType::bondOrder(bond->getBondTypeId()) << endl;
					const auto bondType = bond->getBondTypeId();
					info.order = info.order + BondType::bondOrder(bondType);
					if (bondType == BondType::BondTypeId::AR)
					info.nAromatic++;
				}
				return info;
			};
	orderInfo = reduce<unique_ptr<Bond>, OrderInfo>(bonds, orderInfo, accum);
	double bondOrder = orderInfo.order;
	int nAromatic = orderInfo.nAromatic;

	int doubleBondOrder = round(bondOrder * 2);
	if (AtomType::isCarbonType(getAtomTypeId()) && doubleBondOrder == 9
			&& nAromatic == 3) {
		return 4;
	}
	if (AtomType::isNitrogenType(getAtomTypeId()) && doubleBondOrder == 7
			&& nAromatic == 1) {
		return 3;
	}
	// nitrogen in charged Arginine group w sybyl representation
	if (AtomType::isNitrogenType(getAtomTypeId()) && doubleBondOrder == 9
			&& nAromatic == 3) {
		return 3;
	} else if (AtomType::isSulphurType(getAtomTypeId()) && bondOrder == 3
			&& nAromatic == 2) {
		return 2;
	} else if (AtomType::isOxygenType(getAtomTypeId()) && bondOrder == 3
			&& nAromatic == 2) {
		return 2;
	}
	// Oxygen in carboxylate w sybyl representation
	else if (AtomType::isOxygenType(getAtomTypeId()) && doubleBondOrder == 3
			&& nAromatic == 1) {
		return 1;
	} else if (doubleBondOrder % 2 == 1) {
		auto message = (boost::format("Irregular bond order for atom %s")
				% info()).str();
		throw runtime_error(message);
	}
	return round(bondOrder);
}

const AtomNeighbourhood& Atom::buildNeighbourhood(const Molecule & molecule) {
	neighbourhood = make_unique<AtomNeighbourhood>();
	neighbourhood->create(*this, molecule);

	return *neighbourhood;
}

const int AtomNeighbourhood::countRings() const {
	return rings.size();
}

const bool AtomNeighbourhood::inRingOfSize(size_t size) const {
	for (const auto & ring : rings) {
		if (ring->getAtoms().size() == size) {
			return true;
		}
	}
	return false;
}

const int AtomNeighbourhood::countRingConnections() const {
	int count = 0;
	for (const auto & bond : bonds) {
		if (bond->isInRing()) {
			count++;
		}
	}
	return count;
}

const bool AtomNeighbourhood::isAromatic() const {
	// for the purposes of smarts matching only match ring atoms as aromatic
	return countRings() && countAromaticBonds() > 1;
}

const AtomType::AtomTypeId Atom::checkAtomType(
		const Molecule & molecule) const {

	if (AtomType::isCarbonType(getAtomTypeId())) {
		return checkCarbonType(molecule);
	} else if (AtomType::isNitrogenType(getAtomTypeId())) {
		return checkNitrogenType(molecule);
	} else if (AtomType::isOxygenType(getAtomTypeId())) {
		return checkOxygenType(molecule);
	} else if (AtomType::isSulphurType(getAtomTypeId())) {
		return checkSulphurType(molecule);
	} else if (AtomType::isPhosphorousType(getAtomTypeId())) {
		return checkPhosphorusType(molecule);
	}

	return getAtomTypeId();
}

const AtomTypeId Atom::checkCarbonType(const Molecule & molecule) const {
	assert(AtomType::isCarbonType(getAtomTypeId()));

	auto & neighbours = neighbourhood->getAtomNeighbours();
	auto nNeighbours = static_cast<int>(neighbours.size());
	auto inRing = neighbourhood->getRings().size() > 0;
	auto nSingleNeighbours = neighbourhood->countSingleAndAmideBonds();
	auto nMetalNeighbours = neighbourhood->countMetalNeighbours();
	auto nAromaticNeighbours = neighbourhood->countAromaticBonds();
	auto nDoubleNeighbours = neighbourhood->countDoubleBonds();
	auto nTripleNeighbours = neighbourhood->countTripleBonds();

	if (nNeighbours == 4 && nSingleNeighbours == 4) {
		return AtomTypeId::C3;
	}
	if (nNeighbours == 3 && nSingleNeighbours == 3 && getFormalCharge() == -1) {
		return AtomTypeId::C3;
	}

	if (nNeighbours == 3) {
		// check for arginine carbon atom
		if (nAromaticNeighbours == 3
				|| (nDoubleNeighbours == 1 && nSingleNeighbours == 2)) {
			// lambda will return true if an atom is not nitrogen or
			// does not have three connections
			auto notNpl3Finder =
					[] (const Atom * atom) {
						if (! AtomType::isNitrogenType(atom->getAtomTypeId())) {
							return true;
						}
						if (! (atom->getNeighbourhood().getAtomNeighbours().size()==3)) {
							return true;
						}
						return false;
					};
			if (find_if(neighbours.cbegin(), neighbours.cend(), notNpl3Finder)
					== neighbours.end()) {
				return AtomTypeId::CCAT;
			}
		}
		if (isCarboxylateCarbon(molecule)) {
			return AtomTypeId::C2;
		}

		if (nAromaticNeighbours == 2 && nSingleNeighbours == 1 && inRing
				&& nDoubleNeighbours == 0) {
			return AtomTypeId::CAR;
		}

		// this is unusual!- need to code up Kekulizer also
		if (nAromaticNeighbours == 2 && nDoubleNeighbours == 1 && inRing
				&& nSingleNeighbours == 0) {
			return AtomTypeId::CAR;
		}

		if (nDoubleNeighbours == 1 && nSingleNeighbours == 2) {
			return AtomTypeId::C2;
		}

		if (nAromaticNeighbours == 3) {
			return AtomTypeId::CAR;
		}
	}

	if (nNeighbours == 2 && nTripleNeighbours == 1 && nSingleNeighbours == 1) {
		return AtomTypeId::C1;
	}
	if (nNeighbours == 2 && nDoubleNeighbours == 2) {
		return AtomTypeId::C1;
	}

	// metal ion coordination
	if (nSingleNeighbours == nNeighbours
			&& nNeighbours - nMetalNeighbours == 4) {
		return AtomTypeId::C3;
	}
	if (nSingleNeighbours == nNeighbours - 1 && nDoubleNeighbours == 1
			&& nNeighbours - nMetalNeighbours == 3) {
		return AtomTypeId::C2;
	}
	if (nSingleNeighbours == nNeighbours - 1 && nTripleNeighbours == 1
			&& nNeighbours - nMetalNeighbours == 2) {
		return AtomTypeId::C1;
	}

	if (nNeighbours == nSingleNeighbours && nSingleNeighbours > 4) {
		REPORT(Reporter::WARN) << "Carbon atom with " << nSingleNeighbours
				<< " single bonds";
		return AtomTypeId::C3;
	}

	stringstream ss;
	for (auto bond : neighbourhood->getBonds()) {
		ss << bond->info() << endl;
	}
	REPORT(Reporter::WARN)
			<< "Failed to determine carbon type: bond environment is "
			<< ss.str();
	return AtomTypeId::ATM_NONE;
}

/**
 * Determine if the atom matches the smarts pattern (such that the first term in the
 * smarts matches the atom)
 *
 * @param smarts
 * @param molecule
 * @return
 */
const bool Atom::matchAtom(const string & smarts,
		const Molecule & molecule) const {
	SmartsParser smartsParser;
	auto query = smartsParser.parseSmarts(smarts);
	SubstructureSearch substructureSearch(*query, molecule);
	substructureSearch.fixMapping(0, getAtomNo());
	auto nIsomorphisms = substructureSearch.compare();
	return nIsomorphisms > 0;
}

const bool Atom::isCarboxylateCarbon(const Molecule & molecule) const {
	string smarts("[CX3](~[OX1])(-*)~[OX1]");
	return matchAtom(smarts, molecule);
}

const bool Atom::isPlanarAtom(const Molecule & molecule) const {
	const auto & neighbours = getNeighbourhood().getAtomNeighbours();
	assert(neighbours.size() == 3);
	if (false) {
		// Height above the plane test
		auto atomOop = outOfPlaneHeight(molecule, *this, *neighbours.at(0),
				*neighbours.at(1), *neighbours.at(2));
		return atomOop <= 0.1;
	} else {
		// Use this angle test instead 
		auto midPoint = molecule.getCoord(getAtomNo());
		auto point1 = molecule.getCoord(neighbours.at(0)->getAtomNo());
		auto point2 = molecule.getCoord(neighbours.at(1)->getAtomNo());
		auto point3 = molecule.getCoord(neighbours.at(2)->getAtomNo());
		auto angleSum = .0;
		angleSum += angle(point1, midPoint, point2);
		angleSum += angle(point2, midPoint, point3);
		angleSum += angle(point3, midPoint, point1);
		angleSum *= (180.0 / M_PI);
		return angleSum > 350;
	}
}

const AtomTypeId Atom::checkNitrogenType(const Molecule & molecule) const {
	assert(AtomType::isNitrogenType(getAtomTypeId()));

	const auto & neighbours = neighbourhood->getAtomNeighbours();
	auto nNeighbours = neighbours.size();

	if (nNeighbours == 1) {
		return AtomTypeId::N1;
	}
	if (nNeighbours == 4) {
		return AtomTypeId::N4;
	}

	auto nAromaticNeighbours = neighbourhood->countAromaticBonds();
	auto nSingleNeighbours = neighbourhood->countSingleAndAmideBonds();
	auto nMetalNeighbours = neighbourhood->countMetalNeighbours();
	auto nDoubleNeighbours = neighbourhood->countDoubleBonds();
	auto nTripleNeighbours = neighbourhood->countTripleBonds();

	// common metal ion configuration
	if (nNeighbours - nMetalNeighbours == 4) {
		return AtomTypeId::N4;
	}

	// NPL3 and CCAT- must call checkCarbonType first
	// NPL3 and CCAT
	if (nNeighbours == 3 && nSingleNeighbours == 2) {
		auto ccatFinder = [] (const Atom * atom) {
			return atom->getAtomTypeId() == AtomTypeId::CCAT;
		};
		if (find_if(neighbours.cbegin(), neighbours.cend(), ccatFinder)
				!= neighbours.cend()) {
			return AtomTypeId::NPL3;
		}
	}

	if (nNeighbours == 3 && nDoubleNeighbours == 0) {
		// NAM amide
		if (matchAtom("NC(=O)", molecule)) {
			return AtomTypeId::NAM;
		}

		// NAM or NPL3 for sulfonamide? use NPL3 to be consistent with GOLD, but perhaps NAM is the correct atom type
		if (matchAtom("NS(=O)=O", molecule)) {
			if (molecule.isHasCoordinates()) {
				return isPlanarAtom(molecule) ?
						AtomType::sulfonamideNitrogenType : AtomTypeId::N3;
			} else {
				return AtomType::sulfonamideNitrogenType;
			}
		}

		if (nSingleNeighbours == 3) {

			// NPL3- N single bonded to C or N which is in an aromatic ring or double
			// bonded to a third atom.
			if (matchAtom("N-[#6,#7]=,:*", molecule)) {
				// This is a little more complicated as the N in this environment may be
				// pyrimidal- for 3D molecules also check enviroment to see if N is planar
				// (NPL3) or pyrimial (N3) (e.g. see CCDC multiconformers for 1ywr)
				if (molecule.isHasCoordinates()) {
					return isPlanarAtom(molecule) ?
							AtomTypeId::NPL3 : AtomTypeId::N3;
				} else {
					return AtomTypeId::NPL3;
				}
			} else {
				return AtomTypeId::N3;
			}
		}
	}

	// Nitro group
	if (nNeighbours == 3 && matchAtom("N(~[OX1])~[OX1]", molecule)) {
		return AtomTypeId::N2;
	}

	// Charged Planar Nitrogens
	if (nNeighbours == 3 && nDoubleNeighbours == 1) {
		return AtomTypeId::N2;
	}
	if (nNeighbours == 3 && nAromaticNeighbours == 2) {
		return AtomTypeId::NAR;
	}

	if (nNeighbours == 2 && nDoubleNeighbours == 1) {
		return AtomTypeId::N2;
	}

	if (nNeighbours == 2 && nAromaticNeighbours == 2) {
		return AtomTypeId::NAR;
	}

//	if (nNeighbours == 2 && nSingleNeighbours == 2) {
//		return AtomTypeId::NPL3;
//	}

	// Charged Linear Nitrogen
	if (nNeighbours == 2 && nTripleNeighbours == 1 && nSingleNeighbours == 1) {
		return AtomTypeId::N1;
	}
	if (nNeighbours == 2 && nDoubleNeighbours == 2) {
		return AtomTypeId::N1;
	}

	// this is presumably from a planar nitrogen in a Kekulized ring system.
	// This is taken at face value, but additional checking should still be
	// done. Ring perception etc comes before atom typing
	if (nNeighbours == 3 && nAromaticNeighbours == 3
			&& neighbourhood->countRings() > 0) {
		return AtomTypeId::NAR;
	}

	if (nNeighbours == 5 && nSingleNeighbours == 5 && nMetalNeighbours == 0) {
		// this one is from a smiles in the NCI- don't know what it means!
		REPORT(Reporter::WARN)
				<< "Ntirogen with 5 single bonds (non to metals)!";
		return AtomTypeId::N4;
	}

	if (nNeighbours == 2 && nDoubleNeighbours == 0
			&& matchAtom("NS(=O)=O", molecule)) {
		// this is from the ligand for 1jd0 in the Astex diverse set.
		// The sulfonamide N has only one hydrogen due to metal ion in active site
		REPORT(Reporter::WARN) << "Negatively charged sulfonamide nitrogen!";
		return AtomType::sulfonamideNitrogenType;
	}
	if (nNeighbours == 2 && nDoubleNeighbours == 0
			&& matchAtom("NC(=O)", molecule)) {
		// this is from the ligand in 1n46 in the Astex diverse set.
		// There is a cyclic amide with no hydrogen on the amide N
		REPORT(Reporter::WARN) << "Negatively charged amide nitrogen!";
		return AtomTypeId::NAM;
	}

	stringstream ss;
	for (auto bond : neighbourhood->getBonds()) {
		ss << bond->info() << endl;
	}
	REPORT(Reporter::WARN) << "Failed to determine nitrogen type: for "
			<< info() << " bond environment is " << ss.str();
	return AtomTypeId::ATM_NONE;
}

const AtomTypeId Atom::checkSulphurType(const Molecule & molecule) const {
	assert(AtomType::isSulphurType(getAtomTypeId()));

	auto nOxygens = countExclusivelyBondedO();
	if (nOxygens >= 2) {
		return AtomTypeId::SO2;
	}
	if (nOxygens == 1) {
		return AtomTypeId::SO;
	}

	const auto & neighbours = neighbourhood->getAtomNeighbours();
	auto nNeighbours = neighbours.size();
	if (nNeighbours == 4) {
		return AtomTypeId::S3;
	}

	auto nDoubleNeighbours = neighbourhood->countDoubleBonds();
	auto nSingleNeighbours = neighbourhood->countSingleAndAmideBonds();

	if ((nNeighbours == 1 || nNeighbours == 2 || nNeighbours == 3)
			&& nDoubleNeighbours > 0) {
		return AtomTypeId::S2;
	}
	if (nNeighbours == 2 || nNeighbours == 3) {
		return AtomTypeId::S3;
	}
	if (nNeighbours == 1 && nSingleNeighbours == 1) {
		REPORT(Reporter::WARN) << "Sulphur atom with single connection";
		return AtomTypeId::S3;
	}

	stringstream ss;
	for (auto bond : neighbourhood->getBonds()) {
		ss << bond->info() << endl;
	}
	REPORT(Reporter::WARN)
			<< "Failed to determine sulphur type: bond environment is "
			<< ss.str();
	return AtomTypeId::ATM_NONE;
}

const AtomTypeId Atom::checkPhosphorusType(const Molecule & molecule) const {
	assert(AtomType::isPhosphorousType(getAtomTypeId()));

	return AtomTypeId::P3;
}

const AtomTypeId Atom::checkOxygenType(const Molecule & molecule) const {
	assert(AtomType::isOxygenType(getAtomTypeId()));

	const auto & neighbours = neighbourhood->getAtomNeighbours();
	auto nNeighbours = neighbours.size();
	auto nAromaticNeighbours = neighbourhood->countAromaticBonds();
	auto nSingleNeighbours = neighbourhood->countSingleAndAmideBonds();
	auto nDoubleNeighbours = neighbourhood->countDoubleBonds();
	auto nMetalNeighbours = neighbourhood->countMetalNeighbours();
	auto nTripleNeighbours = neighbourhood->countTripleBonds();

	if (nNeighbours == 2) {
		if (nSingleNeighbours == 2) {
			return AtomTypeId::O3;
		}
		if (neighbourhood->getRings().size() > 0 && nDoubleNeighbours == 1
				&& nSingleNeighbours == 1) {
			return AtomTypeId::O2;
		}
		if (neighbourhood->getRings().size() > 0 && nAromaticNeighbours == 2) {
			return AtomTypeId::OAR;
		}
	}

	if (nNeighbours == 1) {
		// carboxylate
		if (matchAtom("[OX1]~C(-*)~[OX1]", molecule)) {
			return AtomTypeId::OCO2;
		}
		// phosphate
		if (matchAtom("[OX1]~P(-*)~[OX1]", molecule)) {
			return AtomTypeId::OCO2;
		}
		// sulphate
		if (matchAtom("[OX1]:,-S(-*)", molecule)) {
			return AtomTypeId::OCO2;
		}
		if (matchAtom("[OX1]=S(-*)", molecule)) {
			return AtomTypeId::O2;
		}
		if (nDoubleNeighbours == 1) {
			return AtomTypeId::O2;
		}
		// O minus
		if (nSingleNeighbours == 1) {
			return AtomTypeId::O3;
		}
	}

	if (nAromaticNeighbours == 2 && neighbourhood->getRings().size() > 0) {
		REPORT(Reporter::WARN) << "Aromatic oxygen with " << nNeighbours
				<< " connections";
		return AtomTypeId::OAR;
	}

	// common metal ion co-ordinations
	if (nNeighbours == 3 && nSingleNeighbours == 3 && nMetalNeighbours == 1) {
		return AtomTypeId::O3;
	}
	if (nNeighbours == 2 && nDoubleNeighbours == 1 && nMetalNeighbours == 1) {
		return AtomTypeId::O2;
	}

	// O+ triple bond
	if (nNeighbours == 1 && nTripleNeighbours == 1) {
		// unusual case- no SP1 oxygen type
		return AtomTypeId::O2;
	}

	if (nSingleNeighbours > 2) {
		REPORT(Reporter::WARN) << "SP3 Oxygen with " << nNeighbours
				<< " connections";
		return AtomTypeId::O3;
	}
	if (nDoubleNeighbours == 1 && nSingleNeighbours > 0) {
		REPORT(Reporter::WARN) << "SP2 Oxygen with " << nNeighbours
				<< " connections";
		return AtomTypeId::O2;
	}

	stringstream ss;
	for (auto bond : neighbourhood->getBonds()) {
		ss << bond->info() << endl;
	}
	REPORT(Reporter::WARN)
			<< "Failed to determine oxygen type: bond environment is "
			<< ss.str();
	return AtomTypeId::ATM_NONE;
}

const int Atom::countExclusivelyBondedO() const {

	auto accum = [] (int count, const Atom * neighbour) -> int {
		bool ok = AtomType::isOxygenType(neighbour->getAtomTypeId()) &&
		neighbour->getNeighbourhood().getAtomNeighbours().size() == 1;
		auto inc = ok ? 1 : 0;
		return count + inc;
	};
	return reduce<const Atom *, int>(neighbourhood->getAtomNeighbours(), 0,
			accum);
}

} /* namespace GarethMol */
