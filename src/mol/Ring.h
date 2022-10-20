/*
 * Ring.h
 *
 *  Created on: Apr 28, 2014
 *      Author: gjones
 */

#ifndef RING_H_
#define RING_H_

#include <vector>
#include <memory>

#include "Molecule.h"

namespace GarethMol {

class Atom;
using namespace std;

/**
 * Class to represent a ring in a molecular structure
 */
class Ring {

public:

	/**
	 * Create ring from atom indices
	 *
	 * @param molecule
	 * @param atomNos
	 */
	Ring(const Molecule & molecule, const vector<int> & atomNos);
	virtual ~Ring() {
	}

	const std::vector<Atom *>& getAtoms() const {
		return atoms;
	}

	const std::vector<Bond *>& getBonds() const {
		return bonds;
	}

	/**
	 * Return true if this ring contains all sp2 atoms and should thus be planar.
	 *
	 * @return
	 */
	const bool percieveSp2() const;

	const string info() const;

private:

	vector<Atom *> atoms;
	vector<Bond *> bonds;

	Ring(const Ring & other) = delete;
	Ring & operator =(const Ring & other) = delete;

};

} /* namespace GarethMol */

#endif /* RING_H_ */
