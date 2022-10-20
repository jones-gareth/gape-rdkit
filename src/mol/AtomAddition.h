/*
 * AtomAddition.h
 *
 *  Created on: Dec 27, 2015
 *      Author: gareth
 */

#ifndef SRC_MOL_ATOMADDITION_H_
#define SRC_MOL_ATOMADDITION_H_

#include "AtomType.h"
#include "BondType.h"
#include "Molecule.h"
#include "../util/CoordOps.h"
#include <vector>
#include <algorithm>
#include <functional>

namespace GarethMol {

using namespace GarethUtil;
using namespace std;

struct NewAtomInfo {
	const int neighbourAtomNo;
	const AtomType::AtomTypeId & atomTypeId;
	const BondType::BondTypeId & bondTypeId;
	const CoordVector coord;
};

/**
 * This class is used to add atoms and co-ordinates to a molecule
 *
 * A provided adder is used to determine what atoms have new atoms added to them,
 * how many atoms of what type are added.
 *
 * Only terminal atoms can be added.
 */
class AtomAddition {
public:
	AtomAddition(Molecule & mol) :
			molecule(mol) {
	}

	virtual ~AtomAddition() {
	}

	AtomAddition(const AtomAddition & rhs) = delete;
	AtomAddition & operator =(const AtomAddition & rhs) = delete;
	AtomAddition(AtomAddition && rhs) = delete;
	AtomAddition & operator =(AtomAddition && rhs) = delete;

	/**
	 * Adds atoms to the molecule.  Adder tests each atom and returns any
	 * new atoms to be added.
	 *
	 * @param adder
	 */
	void addAtoms(
			function<
					std::vector<NewAtomInfo>(const Atom & existingAtom,
							const Molecule & molecule)> & adder);

	/**
	 * Adds a list of new atom and bond types to a molecule
	 *
	 * @param newAtomsInfo
	 */
	void addAtomsFromList(std::vector<NewAtomInfo> newAtomsInfo);

private:
	Molecule & molecule;
	std::vector<CoordVector> coords;

	void rebuildCoordinateBlock();


};
} /* namespace GarethMol */

#endif /* SRC_MOL_ATOMADDITION_H_ */
