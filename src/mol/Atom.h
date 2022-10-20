/*
 * Atom.h
 *
 *  Created on: Apr 28, 2014
 *      Author: gjones
 */

#ifndef ATOM_H_
#define ATOM_H_

#include <memory>
#include <string>
#include <cassert>

#include "AtomType.h"

namespace GarethMol {

class Bond;
class Molecule;
class Ring;
class AtomNeighbourhood;
class AtomEditKey;
class Angle;
class Torsion;

using namespace std;

class Atom {

public:

	enum class Substructure {
		NONE, RIBOSE, ADENINE, URACIL, BENZENE_, CYTOSINE
	};
	enum class SdfStereo {
		NONE, ODD, EVEN, EITHER
	};

	Atom(int atomNo_, const AtomType & atomType_, const string & name_,
			bool implicit_ = false) :
			atomNo(atomNo_), atomType(&atomType_), name(name_), implicit(
					implicit_) {
		assert(!implicit || atomType->getType() == AtomType::AtomTypeId::H);
	}

	virtual ~Atom() {
	}

	const int getAtomNo() const {
		return atomNo;
	}

	const AtomType & getAtomType() const {
		return *atomType;
	}

	const AtomType::AtomTypeId getAtomTypeId() const {
		return atomType->getType();
	}

	const bool sameAtomType(const Atom & otherAtom) const;

	const string& getName() const {
		return name;
	}

	const SdfStereo getSdfStereo() const {
		return sdfStereo;
	}

	void setSdfStereo(SdfStereo sdfStereo = SdfStereo::NONE) {
		this->sdfStereo = sdfStereo;
	}

	const Substructure getSubstructure() const {
		return substructure;
	}

	void setSubstructure(Substructure substructure = Substructure::NONE) {
		this->substructure = substructure;
	}

	const int getFormalCharge() const {
		return formalCharge;
	}

	void setFormalCharge(int sdfFormalCharge = 0) {
		this->formalCharge = sdfFormalCharge;
	}

	const int getSdfMassDifference() const {
		return sdfMassDifference;
	}

	void setSdfMassDifference(int sdfMassDifference = 0) {
		this->sdfMassDifference = sdfMassDifference;
	}

	const int getSdfValence() const {
		return sdfValence;
	}

	void setSdfValence(int sdfValence = 0) {
		this->sdfValence = sdfValence;
	}

	const bool isFormalChargeSet() const {
		return formalChargeSet;
	}

	void setFormalChargeSet(bool sdfFormalChargeSet = false) {
		this->formalChargeSet = sdfFormalChargeSet;
	}

	const bool isSdfValenceSet() const {
		return sdfValenceSet;
	}

	void setSdfValenceSet(bool sdfValenceSet = false) {
		this->sdfValenceSet = sdfValenceSet;
	}

	/**
	 * determines the current bond order. If the neigbourhood environment has been built,
	 * it is quicker to access from that class
	 *
	 * @return
	 */
	const int totalBondOrder(const Molecule & molecule) const;

	/**
	 * determines the current bond order. If the neigbourhood environment has been built,
	 * it is quicker to access from that class
	 *
	 * @return
	 */
	const int totalBondOrder(const std::vector<unique_ptr<Bond>> & bonds) const;

	const string info() const;

	Atom(const Atom & other) = delete;
	Atom & operator =(const Atom & other) = delete;
	Atom(Atom && other) = delete;
	Atom & operator =(Atom && other) = delete;

	/**
	 * Gets the neighbourhood environment of the atom (assumes it has been built)
	 * @return
	 */
	const AtomNeighbourhood& getNeighbourhood() const {
		return *neighbourhood;
	}

	/**
	 * Predicts the atom type for this atom
	 * @return
	 */
	const AtomType::AtomTypeId checkAtomType(const Molecule & molecule) const;

	/**
	 * builds the local environment of the atom
	 *
	 * @param molecule
	 * @return
	 */
	const AtomNeighbourhood& buildNeighbourhood(const Molecule & molecule);

	const bool isImplicit() const {
		return implicit;
	}

	/**
	 * Allows setting the atom type if you have a key
	 *
	 * @param key
	 * @param atomType
	 */
	void setAtomType(const AtomEditKey & key, const AtomType & atomType) {
		this->atomType = &atomType;
	}

	/**
	 * Returns a non const atom neighbourhood for editing
	 *
	 * @param key
	 * @return
	 */
	AtomNeighbourhood& getNeighbourhoodToEdit(const AtomEditKey & key) {
			return *neighbourhood;
	}

	/**
	 * Matches the molecule to the smarts pattern assuming that the first atom
	 * in the query matches to this atom.
	 *
	 * @param smarts
	 * @param molecule
	 * @return true if the atom and molecule match the pattern
	 */
	const bool matchAtom(const string & smarts,
			const Molecule & molecule) const;

	double getPartialCharge() const {
		return partialCharge;
	}

	void setPartialCharge(double partialCharge = .0) {
		this->partialCharge = partialCharge;
	}

private:
	const int atomNo;
	const AtomType * atomType;
	const string name;
	const bool implicit;

	// information from sd files only
	SdfStereo sdfStereo = SdfStereo::NONE;
	int sdfMassDifference = 0;
	int sdfValence = 0;
	bool sdfValenceSet = false;

	// charges
	int formalCharge = 0;
	double partialCharge = .0;
	bool formalChargeSet = false;

	// neighbour lists
	unique_ptr<AtomNeighbourhood> neighbourhood = nullptr;

	// special substructure
	Substructure substructure = Substructure::NONE;


	// atom type checking functions
	const AtomType::AtomTypeId checkCarbonType(const Molecule & molecule) const;
	const bool isCarboxylateCarbon(const Molecule & molecule) const;
	const AtomType::AtomTypeId checkNitrogenType(
			const Molecule & molecule) const;
	const AtomType::AtomTypeId checkSulphurType(
			const Molecule & molecule) const;
	const AtomType::AtomTypeId checkPhosphorusType(
			const Molecule & molecule) const;
	const AtomType::AtomTypeId checkOxygenType(const Molecule & molecule) const;
	/**
	 * Counts the number of Oxygens bonded to this atom, and to no other atom
	 * @param molecule
	 * @return
	 */
	const int countExclusivelyBondedO() const;

	/**
	 * Returns true if this atom has threee neighbours and is planar
	 */
	const bool isPlanarAtom(const Molecule & molecule) const;

};

/**
 *  A class that contains bond lists and neighbour lists for an atom. Maintains lists of
 *
 *  All neighbours
 *  Atom (real- ie not lone pair or Du) neighbours
 *  Heavy atom neighbours- real and not hydrogen
 *  All bonds
 *  Atom bonds (i.e. bonds to neighbours that are real, not lone pair or Du)
 *  Rings-
 *  I the ring information is not set when the atom neighbourhood is created,
 *  you can use Molecule#updateNeighbourHoodWiithRings once ring perception is done.
 *
 *  The atom neighbourhood will need to be recalculated if
 *  an atom is added or deleted
 *  a real atom is converted to a not-real type or vice versa
 *  a heacy atom is converted to a hydrogen or vice versa
 *  Otherwise atom and bond types may be edited without having to rebuild
 *  the neighbourhood.
 */
class AtomNeighbourhood {

public:
	void create(const Atom & atom, const Molecule & molecule);
	AtomNeighbourhood() {
	}
	virtual ~AtomNeighbourhood() {
	}
	AtomNeighbourhood(const AtomNeighbourhood & rhs) = delete;
	AtomNeighbourhood & operator =(const AtomNeighbourhood & rhs) = delete;
	AtomNeighbourhood(AtomNeighbourhood && rhs) = delete;
	AtomNeighbourhood & operator =(AtomNeighbourhood && rhs) = delete;

	const int countAromaticBonds() const;
	const int countSingleAndAmideBonds() const;
	const int countDoubleBonds() const;
	const int countTripleBonds() const;
	const int totalBondOrder() const;

	const std::vector<const Atom *>& getAllNeighbours() const {
		return allNeighbours;
	}

	const std::vector<const Bond *>& getAtomBonds() const {
		return atomBonds;
	}

	const std::vector<const Atom *>& getAtomNeighbours() const {
		return atomNeighbours;
	}

	const std::vector<const Bond *>& getBonds() const {
		return bonds;
	}

	const std::vector<const Atom *>& getHeavyAtomNeighbours() const {
		return heavyAtomNeighbours;
	}

	const int getNoHydrogens() const {
		return noHydrogens;
	}

	const int countDegree(bool countImplicit = true) const;

	const int countImplicitHydrogens() const;

	const int countRings() const;

	const bool inRingOfSize(size_t size) const;

	const int countRingConnections() const;

	const bool isAromatic() const;

	const int countMetalNeighbours() const;

	void addRing(const Ring * ring) {
		rings.push_back(ring);
	}
	void clearRings() {
		rings.clear();
	}
	const std::vector<const Ring * >& getRings() const {
		return rings;
	}

	const std::vector<const Angle*>& getAngles() const {
		return angles;
	}

	const std::vector<const Torsion*>& getTorsions() const {
		return torsions;
	}

private:
	int noHydrogens = 0;
	std::vector<const Bond *> bonds;
	std::vector<const Bond *> atomBonds;
	std::vector<const Ring *> rings;
	std::vector<const Atom *> allNeighbours;
	std::vector<const Atom *> atomNeighbours;
	std::vector<const Atom *> heavyAtomNeighbours;
	std::vector<const Angle *> angles;
	std::vector<const Torsion *> torsions;

};

/**
 * A key class to allow use of certain editing functions
 */
class AtomEditKey {
	friend class Molecule;
	friend class Huckel;

public:
	virtual ~AtomEditKey() {
	}

	AtomEditKey(const AtomEditKey & rhs) = delete;
	AtomEditKey & operator =(const AtomEditKey & rhs) = delete;
	AtomEditKey(AtomEditKey && rhs) = delete;
	AtomEditKey & operator =(AtomEditKey && rhs) = delete;

private:
	AtomEditKey() {
	}


};
} /* namespace GarethMol */

#endif /* ATOM_H_ */
