/*
 * AtomAddition.cpp
 *
 *  Created on: Dec 27, 2015
 *      Author: gareth
 */

#include "AtomAddition.h"

namespace GarethMol {

using namespace std;

void AtomAddition::addAtoms(
		function<
				std::vector<NewAtomInfo>(const Atom & existingAtom,
						const Molecule & molecule)> & adder) {

	auto atomNo = molecule.nAtoms();
	auto bondNo = molecule.nBonds();
	for (const auto & atom : molecule.getAtoms()) {
		auto newAtomsInfo = adder(*atom, molecule);
		for (const auto & newAtomInfo : newAtomsInfo) {
			const auto & atomType = AtomType::typeFromTypeId(
					newAtomInfo.atomTypeId);
			assert(newAtomInfo.neighbourAtomNo == atom->getAtomNo());
			auto newAtom = make_unique<Atom>(atomNo, atomType,
					atomType.getName());
			const auto & bondType = BondType::typeFromTypeId(
					newAtomInfo.bondTypeId);
			auto bond = make_unique<Bond>(bondNo, bondType, atom.get(), newAtom.get());
			molecule.atoms.push_back(move(newAtom));
			molecule.bonds.push_back(move(bond));
			atomNo++;
			bondNo++;
			coords.push_back(newAtomInfo.coord);
		}
	}

	rebuildCoordinateBlock();
}

void AtomAddition::addAtomsFromList(std::vector<NewAtomInfo> newAtomsInfo) {
	auto atomNo = molecule.nAtoms();
	auto bondNo = molecule.nBonds();
	for (const auto & newAtomInfo : newAtomsInfo) {
		const auto & atomType = AtomType::typeFromTypeId(
				newAtomInfo.atomTypeId);
		const auto & neighbourAtom = molecule.getAtoms().at(
				newAtomInfo.neighbourAtomNo);
		auto newAtom = make_unique<Atom>(atomNo, atomType, atomType.getName());
		const auto & bondType = BondType::typeFromTypeId(
				newAtomInfo.bondTypeId);
		auto bond = make_unique<Bond>(bondNo, bondType, neighbourAtom.get(), newAtom.get());
		molecule.atoms.push_back(move(newAtom));
		molecule.bonds.push_back(move(bond));
		atomNo++;
		bondNo++;
		coords.push_back(newAtomInfo.coord);
	}

	rebuildCoordinateBlock();
}

void AtomAddition::rebuildCoordinateBlock() {
	// rebuild coordinate block
	auto nAtoms = molecule.nAtoms();
	auto newCoords = CoordMatrix(4, nAtoms);
	int start = molecule.coords.cols();
	newCoords.block(0, 0, 4, start) = molecule.coords;
	for (auto coord : coords) {
		newCoords.col(start++) = coord;
	}
	molecule.coords = newCoords;
}

} /* namespace GarethMol */
