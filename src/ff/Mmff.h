/*
 * Mmff.h
 *
 *  Created on: May 15, 2016
 *      Author: gareth
 */

#ifndef SRC_FF_MMFF_H_
#define SRC_FF_MMFF_H_

#include "../mol/Molecule.h"

namespace GarethFF {

using namespace std;
using namespace GarethMol;

class Mmff {
public:
	Mmff(const Molecule & mol) :
			molecule(mol) {
		setupMolecule();
	}

	virtual ~Mmff() {
	}

	Mmff(const Mmff & rhs) = delete;
	Mmff & operator =(const Mmff & rhs) = delete;
	Mmff(Mmff && rhs) = delete;
	Mmff & operator =(Mmff && rhs) = delete;

private:
	const Molecule & molecule;

	void setupMolecule();

};

} /* namespace GarethFF */

#endif /* SRC_FF_MMFF_H_ */
