/*
 * QueryMolecule.h
 *
 *  Created on: Sep 30, 2015
 *      Author: gjones
 */

#ifndef SRC_MOL_QUERYMOLECULE_H_
#define SRC_MOL_QUERYMOLECULE_H_

#include <memory>
#include "Atom.h"
#include "Bond.h"
#include "SmartsAstNode.h"
#include "../util/Array2D.h"

namespace GarethMol {

using namespace GarethUtil;
using namespace std;

class QueryAtom;
class QueryBond;

using QueryAtomPtr = unique_ptr<QueryAtom>;
using QueryBondPtr = unique_ptr<QueryBond>;
using QueryAtomSmartsPtr = unique_ptr<SmartsAstNode<Atom>>;
using QueryBondSmartsPtr = unique_ptr<SmartsAstNode<Bond>>;

/**
 * Class to represent a query atom
 */
class QueryAtom {
public:
	QueryAtom(const int atomNo_, bool matchHeavyAtomOnly_,
			QueryAtomSmartsPtr & smartsAst_) :
			atomNo(atomNo_), matchHeavyAtomOnly(matchHeavyAtomOnly_), smartsAst(
					move(smartsAst_)) {
	}

	virtual ~QueryAtom() {
	}

	QueryAtom(const QueryAtom & rhs) = delete;
	QueryAtom & operator =(const QueryAtom & rhs) = delete;
	QueryAtom(QueryAtom && rhs) = delete;
	QueryAtom & operator =(QueryAtom && rhs) = delete;

	/**
	 * Return true if this atom matches a structure atom
	 *
	 * @param atom
	 * @return
	 */
	bool matchAtom(const Atom & atom, const Molecule & molecule) const;

	const int getAtomNo() const {
		return atomNo;
	}

	const QueryAtomSmartsPtr & getSmartsAst() const {
		return smartsAst;
	}

private:
	const int atomNo;
	const bool matchHeavyAtomOnly;
	/**
	 * AST that is used to match structure atom to query
	 */
	const QueryAtomSmartsPtr smartsAst;
};

/**
 * Class to represent a query bond
 */
class QueryBond {
public:
	QueryBond(const int bondNo_, QueryBondSmartsPtr & smartsAst_,
			const QueryAtom * atom1_, const QueryAtom * atom2_) :
			bondNo(bondNo_), smartsAst(move(smartsAst_)), atom1(atom1_), atom2(
					atom2_) {
	}

	virtual ~QueryBond() {
	}

	QueryBond(const QueryBond & rhs) = delete;
	QueryBond & operator =(const QueryBond & rhs) = delete;
	QueryBond(QueryBond && rhs) = delete;
	QueryBond & operator =(QueryBond && rhs) = delete;

	/**
	 * Matches this query bond to a structure bond
	 * @param bond
	 * @return
	 */
	bool matchBond(const Bond & bond, const Molecule & molecule) const;

	const QueryAtom * getAtom1() const {
		return atom1;
	}

	const QueryAtom * getAtom2() const {
		return atom2;
	}

	const int getBondNo() const {
		return bondNo;
	}

private:
	const int bondNo;
	/**
	 * AST that is used to match the structure bond to the query
	 */
	const QueryBondSmartsPtr smartsAst;
	const QueryAtom * atom1, *atom2;
};

/**
 * A class to represent a query molecule.
 *
 */
class QueryMolecule {
public:
	QueryMolecule(const string & name_, vector<QueryAtomPtr> & atoms_,
			 vector<QueryBondPtr> & bonds_);

	virtual ~QueryMolecule() {
	}

	QueryMolecule(const QueryMolecule & other) = delete;
	QueryMolecule & operator =(const QueryMolecule & other) = delete;
	QueryMolecule(Atom && other) = delete;
	QueryMolecule & operator =(QueryMolecule && other) = delete;

	const vector<QueryAtomPtr>& getAtoms() const {
		return atoms;
	}

	const vector<QueryBondPtr>& getBonds() const {
		return bonds;
	}

	const string& getName() const {
		return name;
	}

	const QueryBond * getBond(const int atom1, const int atom2) const;

private:
	const string name;
	const vector<QueryAtomPtr> atoms;
	const vector<QueryBondPtr> bonds;
	Array2D<QueryBond *> bondTable;

};

} /* namespace GarethMol */

#endif /* SRC_MOL_QUERYMOLECULE_H_ */
