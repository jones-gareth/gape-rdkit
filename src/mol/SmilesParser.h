/*
 * Classes for parsing smiles and converting them to molecules
 *
 * Does not parse the following:cat
 * Isotropic specification
 * Bond stereochemistry
 * Atom chirality
 *
 * SmilesParser.h
 *
 *  Created on: Sep 22, 2015
 *      Author: gjones
 */

#ifndef SRC_MOL_SMILESPARSER_H_
#define SRC_MOL_SMILESPARSER_H_

#include <set>
#include <string>
#include <map>
#include <memory>
#include "Molecule.h"

namespace GarethMol {

/**
 * Hold a ring closures in a smiles
 */
class RingClosureInformation {
public:
	RingClosureInformation(const int closureNo_, const Atom * atom_,
			const BondType::BondTypeId bondType_) :
			closureNo(closureNo_), atom(atom_), bondType(bondType_) {
	}

	virtual ~RingClosureInformation() {
	}

	RingClosureInformation(const RingClosureInformation & rhs) = delete;
	RingClosureInformation & operator =(const RingClosureInformation & rhs) = delete;
	RingClosureInformation(RingClosureInformation && rhs) = delete;
	RingClosureInformation & operator =(RingClosureInformation && rhs) = delete;

	const Atom * getAtom() const {
		return atom;
	}

	const BondType::BondTypeId getBondType() const {
		return bondType;
	}

	const int getClosureNo() const {
		return closureNo;
	}

private:

	const int closureNo;
	const Atom * atom;
	const BondType::BondTypeId bondType;
};

/**
 * Main parsing class
 */
class SmilesParser {
public:
	SmilesParser();
	virtual ~SmilesParser();

	SmilesParser(const SmilesParser & rhs) = delete;
	SmilesParser & operator =(const SmilesParser & rhs) = delete;
	SmilesParser(SmilesParser && rhs) = delete;
	SmilesParser & operator =(SmilesParser && rhs) = delete;

	/**
	 * Main routine to convert a smiles to a Molecule.
	 *
	 * @param smiles_
	 * @return
	 */
	unique_ptr<Molecule> parseSmiles(const string & smiles_);

private:
	// our current parsing position
	size_t position = 0;
	// input smiles
	string smiles { "" };
	// Atoms and bonds
	vector<unique_ptr<Atom>> atoms { };
	vector<unique_ptr<Bond>> bonds { };
	bool firstAtomRead = false;
	Atom * currentAtom = nullptr;
	size_t atomNo = 0, bondNo = 0;
	BondType::BondTypeId currentBond = BondType::BondTypeId::UNK;
	vector<Atom *> openBranches { };
	map<size_t, unique_ptr<RingClosureInformation>> ringOpenings { };
	set<size_t> aromaticAtomNumbers { };
	bool currentAtomAromatic = false;

	/**
	 * Parse an atom type
	 *
	 * @param organicSubset
	 * @return
	 */
	const AtomType & parseAtomType(const bool organicSubset);

	/**
	 *
	 * @param pos
	 * @return the character from the smiles string at the position
	 */
	char getCharacter(size_t pos) const;

	/**
	 *
	 * @return the current characterin the smiles
	 */
	char getCurrentCharacter() const;

	/**
	 * Throw a parsing error
	 * @param message
	 */
	void parseError(const string & message);

	/**
	 * Format a parse message with the current position
	 *
	 * @param message
	 * @return
	 */
	const string parseMessage(const string & message);

	/**
	 * Handle an atom entry
	 */
	void parseAtomEntry();

	/**
	 * Handle a branch opening
	 */
	void parseOpenBranch();

	/**
	 * Handle a bond
	 */
	void parseBond();

	/**
	 * Parse an atom definition
	 */
	void parseAtom();

	/**
	 * Parse a ring closure
	 *
	 * @return
	 */
	bool parseRingClosure();

	/**
	 * Parse atom attributes (charge etc)
	 */
	void parseAtomAttributes();

	/**
	 * Handle a brance closure
	 */
	void parseCloseBranch();

	/**
	 * Add a number of hydrogens to the current atom
	 *
	 * @param noHydrogens
	 */
	void addHydrogens(int noHydrogens);

	/**
	 * Add hydrogens to the atom
	 *
	 * @param atom
	 * @param noHydrogens
	 * @param implicit
	 */
	void addHydrogens(const Atom * atom, int noHydrogens, const bool implicit);

	/**
	 * Fill valence adding hydrogens
	 */
	void addHydrogens();

	/**
	 * Add an atom to the structure under construction
	 *
	 * @param type
	 * @param implicit  set true if this is an implicit hydrogen
	 * @return
	 */
	Atom * addAtom(const AtomType & type, const bool implicit = false);

	/**
	 * Adds a bond to the structure under construction
	 *
	 * @param typeId
	 * @param atom1
	 * @param atom2
	 * @return
	 */
	Bond * addBond(const BondType::BondTypeId typeId, const Atom * atom1,
			const Atom * atom2);

	/**
	 *
	 * @param atom
	 * @return true if this atom is identified as aromatic in the smiles string
	 */
	const bool isAromatic(const Atom * atom) const;

	/**
	 * Gets an integer from the current position in the smiles string
	 * @return
	 */
	const int getInteger();

	/**
	 * Determines the number of times the test character is at the current position
	 * in the smiles string
	 *
	 * @param testChar
	 * @return
	 */
	const int countRepeatingCharacter(const char testChar);

	/**
	 * Get any charge definition from the current position in the smiles string
	 * @param chargeChar
	 * @return
	 */
	const int getCharge(const char chargeChar);
};

} /* namespace GarethMol */

#endif /* SRC_MOL_SMILESPARSER_H_ */
