/*
 * MoleculeSssr.cpp
 *
 *  Created on: Nov 3, 2015
 *      Author: gjones
 */

#include "MoleculeSssr.h"
#include "Ring.h"
#include "../util/Util.h"
#include "../util/Reporter.h"

using namespace GarethUtil;

namespace GarethMol {

void MoleculeSssr::doSssr() {
	Sssr<MoleculeSssr> sssr(*this);
	molecule.createBondTableIfAbsent();
	vector<vector<int>> ringsVector = sssr.doSssr();
	molecule.clearRings();

	for (const auto & ringAtomIds : ringsVector) {
		REPORT(Reporter::DEBUG) << "Adding ring of size " << ringAtomIds.size();
		// create ring
		auto ring = make_unique<Ring>(molecule, ringAtomIds);
		// add ring to molecule
		molecule.addRing(ring);
	}

}

const vector<int> MoleculeSssr::getNodeNeighbours(const int atomNo) const {
	auto filterFunction = [atomNo] (const Bond & bond) {
		return bond.getAtom1().getAtomNo() == atomNo ||
		bond.getAtom2().getAtomNo() == atomNo;
	};
	auto atomBonds = filterUniquePtrListToRawPtrList<Bond>(molecule.getBonds(),
			filterFunction);
	auto mapFunction = [atomNo] (  Bond * bond) -> int{
		if (bond->getAtom1().getAtomNo() == atomNo) {
			return bond->getAtom2().getAtomNo();
		}
		if (bond->getAtom2().getAtomNo() == atomNo) {
			return bond->getAtom1().getAtomNo();
		}
		assert(false);
	};
	return mapToNewList<Bond *, int>(atomBonds, mapFunction);
}

} /* namespace GarethMol */
