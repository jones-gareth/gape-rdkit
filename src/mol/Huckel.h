/*
 * Huckel.h
 *
 *  Created on: Jan 21, 2016
 *      Author: gjones
 */

#ifndef SRC_MOL_HUCKEL_H_
#define SRC_MOL_HUCKEL_H_


#include "Atom.h"
#include "Molecule.h"
#include <set>

/**
 * Performs aromatic ring perception using Huckels rule.
 *
 * The molecule has to have SSSR determinied and atom neighbourhoods added.
 * Also, types assigned if loaded from an SDF file.
 */
namespace GarethMol {

using RingSystem = std::set<int>;

class Huckel {
public:
	Huckel(const Molecule & mol): molecule(mol) {
	}

	virtual ~Huckel() {
	}

	Huckel(const Huckel & rhs) = delete;
	Huckel & operator =(const Huckel & rhs) = delete;
	Huckel(Huckel && rhs) = delete;
	Huckel & operator =(Huckel && rhs) = delete;

	void findAromaticRings();
	const string info() const;

private:
	const Molecule & molecule;
	std::vector<RingSystem> ringSystems{};

	void findRingSystems();
	bool isAromaticRingSystem(const RingSystem &ringSystem) const;
	bool sp2ok(const Atom & atom, const RingSystem &ringSystem) const;
	void setRingAromatic(const Ring & ring);
};

} /* namespace GarethMol */

#endif /* SRC_MOL_HUCKEL_H_ */
