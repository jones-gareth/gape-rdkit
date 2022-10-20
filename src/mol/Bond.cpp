/*
 * Bond.cpp
 *
 *  Created on: Apr 28, 2014
 *      Author: gjones
 */

#include "Bond.h"
#include "Molecule.h"

namespace GarethMol {

using BondTypeId = BondType::BondTypeId;
using AtomTypeId = AtomType::AtomTypeId;

bool Bond::sameBondType(const Bond & otherBond) const {
	return getBondType().getType() == otherBond.getBondType().getType();
}

const string Bond::info() const {
	return to_string(bondNo + 1) + ":" + getBondType().getName();
}

/**
 * Return true if atom1 is an amide nitrogen and atom2 is the amide carbon or sulphur
 *
 * @param atom1
 * @param atom2
 * @param molecule
 * @return
 */
static const bool amideBond(const Atom & atom1, const Atom & atom2,
		const Molecule & molecule) {
	if (atom1.getAtomTypeId() == AtomTypeId::NAM
			&& atom2.getAtomTypeId() == AtomTypeId::C2
			&& atom2.matchAtom("C(=O)N", molecule)) {
		return true;
	}

	// this currently shouldn't get called as sulphonamide nitrogen atoms are typed as
	// NPL3 (unless the N is common to an amide and sulfonamide group
	if (atom1.getAtomTypeId() == AtomTypeId::NAM
			&& AtomType::isSulphurType(atom2.getAtomTypeId())
			&& atom2.matchAtom("S(=O)(=O)N", molecule)) {
	    if(AtomType::sulfonamideNitrogenType == AtomType::AtomTypeId::NAM) {
	        assert(atom1.matchAtom("NC(=O)", molecule));
	        return false;
	    }
		return true;
	}

	return false;

}

const BondTypeId Bond::sybylBondType(const Molecule & molecule) const {
		if (atom1->getAtomTypeId() == AtomTypeId::OCO2
				&& AtomType::isRealAtom(atom2->getAtomTypeId())) {
			return BondTypeId::AR;
		}
		if (atom2->getAtomTypeId() == AtomTypeId::OCO2
				&& AtomType::isRealAtom(atom1->getAtomTypeId())) {
			return BondTypeId::AR;
		}
		if (atom1->getAtomTypeId() == AtomTypeId::CCAT
				&& AtomType::isNitrogenType(atom2->getAtomTypeId())) {
			return BondTypeId::AR;
		}
		if (atom2->getAtomTypeId() == AtomTypeId::CCAT
				&& AtomType::isNitrogenType(atom1->getAtomTypeId())) {
			return BondTypeId::AR;
		}

		return bondType->getType();
}

const BondTypeId Bond::checkBondType(const Molecule & molecule) const {

	if (false) {
		// using the sybybl bond types screws up formal charge assignement and
		// SD file handling- a separate syblyBondType method is provided to convert
		//to sybyl bond types.  Do we need to have something to deal with when
		// reading MOL2 files?
		if (atom1->getAtomTypeId() == AtomTypeId::OCO2
				&& AtomType::isRealAtom(atom2->getAtomTypeId())) {
			return BondTypeId::AR;
		}
		if (atom2->getAtomTypeId() == AtomTypeId::OCO2
				&& AtomType::isRealAtom(atom1->getAtomTypeId())) {
			return BondTypeId::AR;
		}
		if (atom1->getAtomTypeId() == AtomTypeId::CCAT
				&& AtomType::isNitrogenType(atom2->getAtomTypeId())) {
			return BondTypeId::AR;
		}
		if (atom2->getAtomTypeId() == AtomTypeId::CCAT
				&& AtomType::isNitrogenType(atom1->getAtomTypeId())) {
			return BondTypeId::AR;
		}
	}
	if (amideBond(*atom1, *atom2, molecule)) {
		return BondTypeId::AM;
	}
	if (amideBond(*atom2, *atom1, molecule)) {
		return BondTypeId::AM;
	}

	return getBondType().getType();
}

}
/* namespace GarethMol */
