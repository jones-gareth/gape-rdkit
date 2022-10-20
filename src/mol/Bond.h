/*
 * Bond
 *
 *  Created on: Apr 28, 2014
 *      Author: gjones
 */

#ifndef BOND_H_
#define BOND_H_

#include <memory>

#include "BondType.h"

namespace GarethMol {

class Ring;
class Atom;
class Molecule;
class BondEditKey;

using namespace std;

/**
 * Class to represent a chemical bond
 */
class Bond {

public:

	/**
	 * Enumeration for sbond stereo settings in an SDF file
	 */
	enum class SdfStereo {
		NONE, UP, DOWN, EITHER, CIS_TRANS
	};

	/**
	 * Enumeration for rotatable bonds
	 */
	enum class Rotatable {
		NO, FULL, FLIP
	};

	/**
	 * Construct a bond between two atoms.
	 *
	 * @param bondNo_ bond number
	 * @param bt bondf type
	 * @param a1 atom 1
	 * @param a2 atom 2
	 */
	explicit Bond(const int bondNo_, const BondType & bt, const Atom * a1,
			const Atom * a2) :
			bondNo(bondNo_), bondType(&bt), atom1(a1), atom2(a2) {
	}

	virtual ~Bond() {
	}

	// getter and setters

	bool isOutput() const {
		return output;
	}

	void setOutput(bool output = true) {
		this->output = output;
	}

	SdfStereo getStereo() const {
		return stereo;
	}

	void setStereo(const SdfStereo stereo_ = SdfStereo::NONE) {
		stereo = stereo_;
	}

	const Rotatable & getRotatable() const {
		return rotatable;
	}

	void setRotatable(const Rotatable rotatable_ = Rotatable::NO) {
		rotatable = rotatable_;
	}

	const int getBondNo() const {
		return bondNo;
	}

	const BondType & getBondType() const {
		return *bondType;
	}

	const BondType::BondTypeId getBondTypeId() const {
		return bondType->getType();
	}

	const BondType::BondTypeId getBondTypeId(const BondType & bondType) const {
		return bondType.getType();
	}

	const Atom & getAtom1() const {
		return *atom1;
	}

	const Atom & getAtom2() const {
		return *atom2;
	}

	const Atom * getAtom1Ptr() const {
		return atom1;
	}

	const Atom * getAtom2Ptr() const {
		return atom2;
	}

	bool sameBondType(const Bond & otherBond) const;

	bool isInRing() const {
		return inRing;
	}

	void setInRing(bool inRing = false) {
		this->inRing = inRing;
	}

	/**
	 * Checks the bond type
	 *
	 * @param molecule
	 * @return
	 */
	const BondType::BondTypeId checkBondType(const Molecule & molecule) const;

	/**
	 * Returns the Tripos bond type.  This is "aromatic" for guanidinium
	 * and carboxylate
	 *
	 * @param molecule
	 * @return
	 */
	const BondType::BondTypeId sybylBondType(const Molecule & molecule) const;

	/**
	 * Return an informational string.
	 *
	 * @return
	 */
	const string info() const;

	void setBondType(const BondEditKey & key, const BondType & bondType) {
		this->bondType = &bondType;
	}

private:

	const int bondNo;
	const BondType * bondType;
	const Atom * atom1;
	const Atom * atom2;
	bool inRing = false;

	std::vector<Ring *> rings{};

	SdfStereo stereo { SdfStereo::NONE };
	Rotatable rotatable { Rotatable::NO };
	bool output = true;

	Bond(const Bond & other) = delete;
	Bond & operator =(const Bond & other) = delete;

};

/**
 * A key class to allow editing of bond types
 */
class BondEditKey {
	friend class Molecule;
	friend class Huckel;

public:

	virtual ~BondEditKey() {
	}

	BondEditKey(const BondEditKey & rhs) = delete;
	BondEditKey & operator =(const BondEditKey & rhs) = delete;
	BondEditKey(BondEditKey && rhs) = delete;
	BondEditKey & operator =(BondEditKey && rhs) = delete;

private:
	BondEditKey() {
	}
};
} /* namespace GarethMol */

#endif /* BOND_H_ */
