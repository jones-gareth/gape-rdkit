/*
 * QueryMolecule.cpp
 *
 *  Created on: Sep 30, 2015
 *      Author: gjones
 */

#include "QueryMolecule.h"

namespace GarethMol {

QueryMolecule::QueryMolecule(const string & name_, vector<QueryAtomPtr> & atoms_,
		 vector<QueryBondPtr> & bonds_) :
		name(name_), atoms(move(atoms_)), bonds(move(bonds_)), bondTable(
				atoms.size(), atoms.size()) {
	for (auto & bond : bonds) {
		bondTable.get(bond->getAtom1()->getAtomNo(),
				bond->getAtom2()->getAtomNo()) = bond.get();
		bondTable.get(bond->getAtom2()->getAtomNo(),
				bond->getAtom1()->getAtomNo()) = bond.get();
	}
}

bool QueryAtom::matchAtom(const Atom & atom, const Molecule & molecule) const {
	if (matchHeavyAtomOnly && atom.getAtomTypeId() == AtomType::AtomTypeId::H)
		return false;
	return smartsAst->evaluate(atom, molecule);
}

bool QueryBond::matchBond(const Bond & bond, const Molecule & molecule) const {
	return smartsAst->evaluate(bond, molecule);
}

const QueryBond * QueryMolecule::getBond(const int atom1,
		const int atom2) const {
	auto bond = bondTable.get(atom1, atom2);
	return bond;
}

} /* namespace GarethMol */
