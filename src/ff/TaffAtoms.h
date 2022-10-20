/*
 * TaffAtoms.h
 *
 *  Created on: Mar 8, 2016
 *      Author: gareth
 */

#ifndef SRC_FF_TAFFATOMS_H_
#define SRC_FF_TAFFATOMS_H_

#include <map>
#include <memory>
#include <vector>

#include "../mol/AtomType.h"

namespace GarethMol {
class Molecule;
} /* namespace GarethMol */

namespace GarethFF {

using namespace GarethMol;
using namespace std;

class TaffAtom {

	using AtomTypeId = GarethMol::AtomType::AtomTypeId;

public:

	TaffAtom(const double r_, const double k_, const AtomTypeId id_) :
			r(r_), k(k_), atomTypeId(id_) {
	}

	virtual ~TaffAtom() {
	}

	TaffAtom(const TaffAtom & rhs) = delete;
	TaffAtom & operator =(const TaffAtom & rhs) = delete;
	TaffAtom(TaffAtom && rhs) = delete;
	TaffAtom & operator =(TaffAtom && rhs) = delete;

	const AtomTypeId getAtomTypeId() const {
		return atomTypeId;
	}

	const double getK() const {
		return k;
	}

	const double getR() const {
		return r;
	}

private:
	const double r;
	const double k;
	const AtomTypeId atomTypeId;


};

class TaffAtoms {
public:

	static const TaffAtoms & getInstance();

	virtual ~TaffAtoms() {
	}

	TaffAtoms(const TaffAtoms & rhs) = delete;
	TaffAtoms & operator =(const TaffAtoms & rhs) = delete;
	TaffAtoms(TaffAtoms && rhs) = delete;
	TaffAtoms & operator =(TaffAtoms && rhs) = delete;

	map<int, TaffAtom *> parametersForMolecule(const Molecule & mol) const;
private:


	TaffAtoms() {
	}

	static std::vector<unique_ptr<TaffAtom>> & taffAtomDefinitions() ;
};

} /* namespace GarethFF */

#endif /* SRC_FF_TAFFATOMS_H_ */
