/*
 * Huckel.cpp
 *
 *  Created on: Jan 21, 2016
 *      Author: gjones
 */

#include "Huckel.h"

#include <stddef.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <memory>
#include <utility>
#include <vector>

#include "Ring.h"
#include "../util/Reporter.h"
#include "../util/Util.h"

namespace GarethMol {

using namespace std;
using namespace GarethUtil;

/**
 * Converts a Ring to a RingSystem set of atom numbers
 *
 * @param ring
 * @return
 */
static RingSystem toRingSystem(const Ring & ring) {
	RingSystem rs;
	for (const auto & atom : ring.getAtoms()) {
		rs.insert(atom->getAtomNo());
	}
	return rs;
}

/**
 * Returns true if ring1 and ring2 share any atoms
 *
 * @param ring1
 * @param ring2
 * @return
 */
static bool sameRingSystem(const RingSystem & ring1, const RingSystem & ring2) {
	for (auto no1 : ring1) {
		if (ring2.count(no1)) {
			return true;
		}
	}
	return false;
}

/**
 * Merges two ring systems and returns the union
 *
 * @param ring1
 * @param ring2
 * @return
 */
static RingSystem mergeRings(const RingSystem & ring1,
		const RingSystem & ring2) {
	RingSystem merged;
	set_union(ring1.cbegin(), ring1.cend(), ring2.cbegin(), ring2.cend(),
			insert_iterator<set<int>>(merged, merged.end()));
	return merged;
}

/**
 * Return true if the ring is wholly contained in the ring system
 *
 * @param ring
 * @param ringSystem
 * @return
 */
static bool inRingSystem(const Ring & ring, const RingSystem &ringSystem) {
	for (const auto & atom : ring.getAtoms()) {
		if (!ringSystem.count(atom->getAtomNo())) {
			return false;
		}
	}
	return true;
}

void Huckel::findRingSystems() {

	std::vector<RingSystem> rings;
	for (const auto & ring : molecule.getRings()) {
		if (ring->percieveSp2()) {
			rings.push_back(toRingSystem(*ring));
		}
	}

	while (true) {

		auto nRings = rings.size();
		size_t ring1, ring2;
		bool foundRings = false;
		for (auto i = 0ul; i < nRings; i++) {
			for (auto j = i + 1; j < nRings; j++) {
				if (sameRingSystem(rings.at(i), rings.at(j))) {
					ring1 = i;
					ring2 = j;
					foundRings = true;
					goto FOUND;
				}
			}
		}
		FOUND: if (!foundRings) {
			ringSystems = rings;

			REPORT(Reporter::DEBUG) << "Found ring systems " << info();
			return;
		}
		std::vector<RingSystem> newRings;
		newRings.push_back(mergeRings(rings.at(ring1), rings.at(ring2)));

		for (auto i = 0ul; i < nRings; i++) {
			if (i != ring1 && i != ring2) {
				newRings.push_back(rings.at(i));
			}
		}
		rings = newRings;
	}
}

const string Huckel::info() const {
	ostringstream ss;
	ss << "[";
	for (const auto & ring : ringSystems) {
		std::vector<int> atoms(ring.cbegin(), ring.cend());
		auto ids = mapToNewList<int>(atoms,
				[](const int atomNo) {return atomNo+1;});
		ss << "[";
		ss << collectionToString(ids, ",");
		ss << "]";
	}
	ss << "]";
	return ss.str();
}

void Huckel::setRingAromatic(const Ring & ring) {
	using AtomTypeId = AtomType::AtomTypeId;
	AtomEditKey atomEditKey;
	for (auto & atom : ring.getAtoms()) {
		AtomTypeId newType = AtomTypeId::DU;
		const AtomTypeId & currentType = atom->getAtomTypeId();
		if (AtomType::isOxygenType(currentType)) {
			newType = AtomTypeId::OAR;
		} else if (AtomType::isCarbonType(currentType)) {
			newType = AtomTypeId::CAR;
		} else if (AtomType::isNitrogenType(currentType)) {
			newType = AtomTypeId::NAR;
		}

		if (newType != AtomTypeId::DU && newType != currentType) {
			const auto & newAtomType = AtomType::typeFromTypeId(newType);
			auto & editAtom = molecule.getAtom(atom->getAtomNo());
			REPORT(Reporter::DEBUG) << "Set atom type for atom "
					<< atom->info() << " to type " << newAtomType.getName();
			editAtom.setAtomType(atomEditKey, newAtomType);
		}
	}

	const auto &newBondType = BondType::typeFromTypeId(
			BondType::BondTypeId::AR);
	BondEditKey bondEditKey;
	for (auto & bond : ring.getBonds()) {
		auto & editBond = molecule.getBond(bond->getBondNo());
		REPORT(Reporter::DEBUG) << "Set bond type for bond " << bond->info()
				<< " to type " << newBondType.getName();
		if (bond->getBondTypeId() != BondType::BondTypeId::AR) {
			editBond.setBondType(bondEditKey, newBondType);
		}
	}
}

void Huckel::findAromaticRings() {
	findRingSystems();

	for (const auto & ringSystem : ringSystems) {
		if (isAromaticRingSystem(ringSystem)) {

			for (const auto & ring : molecule.getRings()) {
				if (inRingSystem(*ring, ringSystem)) {
					REPORT(Reporter::DEBUG) << " ring " << ring->info()
							<< " is aromatic";
					setRingAromatic(*ring);
				}
			}
		}
	}
}

bool Huckel::isAromaticRingSystem(const RingSystem &ringSystem) const {
	REPORT(Reporter::DEBUG) << "Checking ring system for aromaticity";
	auto nElectrons = 0;

	// determine number of electons that can contribute to aromatic system
	for (const auto atomNo : ringSystem) {
		const auto & atom = molecule.getAtom(atomNo);
		auto typeId = atom.getAtomTypeId();

		if (AtomType::isPhosphorousType(typeId)) {
			REPORT(Reporter::DEBUG) << "Atom " << atom.info()
					<< " is not aromatic";
			return false;
		}

		else if (AtomType::isOxygenType(typeId)) {
			nElectrons += 2;
		}

		else if (AtomType::isSulphurType(typeId)) {
			if (sp2ok(atom, ringSystem)) {
				nElectrons += 2;
			} else {
				REPORT(Reporter::DEBUG) << "Atom " << atom.info()
						<< " is not aromatic";
				return false;
			}
		}

		else if (AtomType::isCarbonType(typeId)) {
			if (sp2ok(atom, ringSystem)) {
				nElectrons++;
			} else if (atom.getNeighbourhood().getAtomNeighbours().size() == 2
					&& atom.getNeighbourhood().countSingleAndAmideBonds()
							== 2) {
				// -1 carbon, with two single bonds can contribute 1 electron
				nElectrons++;
			} else {
				REPORT(Reporter::DEBUG) << "Atom " << atom.info()
						<< " is not aromatic";
				return false;
			}
		}

		else if (AtomType::isNitrogenType(typeId)) {
			if (sp2ok(atom, ringSystem)) {
				nElectrons++;
			} else if (atom.getNeighbourhood().getAtomNeighbours().size() == 3
					&& atom.getNeighbourhood().countSingleAndAmideBonds()
							== 3) {
				nElectrons += 2;
			} else {
				REPORT(Reporter::DEBUG) << "Atom " << atom.info()
						<< " is not aromatic";
				return false;
			}
		}

	}

	// 4n +2 rule
	int test = nElectrons - 2;
	if (test % 4 == 0) {
		REPORT(Reporter::DEBUG) << "Ring system passes Huckels rule";
		return true;
	}

	REPORT(Reporter::DEBUG) << "Ring system fails Huckels test";
	return false;
}

bool Huckel::sp2ok(const Atom & atom, const RingSystem &ringSystem) const {

	struct BondCounts {
		int nSingle;
		int nAromatic;
		int nDouble;
		int nTotal;
	};

	BondCounts bondCounts = { 0, 0, 0 };

	auto accum =
			[& atom, & ringSystem](BondCounts & info, const Bond * bond) {
				int atomNo = atom.getAtomNo();
				int otherAtomNo = bond->getAtom1().getAtomNo()==atomNo ? bond->getAtom2().getAtomNo() : bond->getAtom1().getAtomNo();
				assert (otherAtomNo != atomNo);
				if (ringSystem.count(otherAtomNo)>0) {
					const auto bondType = bond->getBondTypeId();
					if (bondType == BondType::BondTypeId::SINGLE) {
						info.nSingle++;
					}
					if (bondType == BondType::BondTypeId::DOUBLE) {
						info.nDouble++;
					}
					if (bondType == BondType::BondTypeId::AR) {
						info.nAromatic++;
					}
					info.nTotal++;
				}
				return info;
			};

	bondCounts = reduce<const Bond *, BondCounts>(
			atom.getNeighbourhood().getAtomBonds(), bondCounts, accum);

	if (bondCounts.nTotal == 2) {
		if (bondCounts.nSingle == 1 && bondCounts.nDouble == 1)
			return true;
		if (bondCounts.nAromatic == 2)
			return true;
	} else if (bondCounts.nTotal == 3) {
		if (bondCounts.nSingle == 2 && bondCounts.nDouble == 1)
			return true;
		if (bondCounts.nAromatic == 3)
			return true;

	}
	REPORT(Reporter::DEBUG) << "Atom " << atom.info()
			<< "is not aromatic in current ring system";
	return false;
}

} /* namespace GarethMol */
