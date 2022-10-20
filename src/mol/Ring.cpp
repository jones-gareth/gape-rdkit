/*
 * Ring.cpp
 *
 *  Created on: Apr 28, 2014
 *      Author: gjones
 */

#include "Ring.h"
#include "../util/Reporter.h"
#include "../util/Util.h"

namespace GarethMol {

using namespace std;
using namespace GarethUtil;

Ring::Ring(const Molecule & molecule, const vector<int> & atomNos) {
	assert(atoms.size() == 0);
	assert(bonds.size() == 0);
	atoms.reserve(atomNos.size());
	bonds.reserve(atomNos.size());

	auto lastAtomNo = *(atomNos.cend() - 1);
	for (auto iter = atomNos.cbegin(); iter < atomNos.cend(); ++iter) {
		auto atomNo = *iter;
		auto & atom = molecule.getAtoms().at(atomNo);
		atoms.push_back(atom.get());
		auto bond = molecule.getBond(lastAtomNo, atomNo);
		assert(bond != nullptr);
		bond->setInRing(true);
		bonds.push_back(bond);
		lastAtomNo = atomNo;
	}
}

const bool Ring::percieveSp2() const {

	for (const auto & atom : atoms) {
		const auto & type = atom->getAtomType();
		if (type.getGeometry() != AtomType::Geometry::TRI) {

			if (type.getType() == AtomType::AtomTypeId::S3) {
				continue;
			}

			if (AtomType::isOxygenType(type.getType())) {
				continue;
			}

			REPORT(Reporter::DEBUG) << "Atom " << atom->info()
					<< " is not SP3";
			return false;
		}
	}

	return true;
}

const string Ring::info() const {
	using AtomPtr = Atom *;
	auto ids = mapToNewList<Atom *, int>(atoms, [] (const AtomPtr & atom) {
		return atom->getAtomNo() + 1;
	});
	return "["+collectionToString(ids, ",")+"]";
}

} /* namespace GarethMol */
