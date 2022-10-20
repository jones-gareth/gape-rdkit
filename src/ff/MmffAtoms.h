/*
 * MmffAtoms.h
 *
 *  Created on: May 15, 2016
 *      Author: gareth
 */

#ifndef SRC_FF_MMFFATOMS_H_
#define SRC_FF_MMFFATOMS_H_

#include <memory>
#include <vector>

namespace GarethFF {

using namespace std;

class MmffAtom {
public:
	MmffAtom() {
	}

	virtual ~MmffAtom() {
	}

	MmffAtom(const MmffAtom & rhs) = delete;
	MmffAtom & operator =(const MmffAtom & rhs) = delete;
	MmffAtom(MmffAtom && rhs) = delete;
	MmffAtom & operator =(MmffAtom && rhs) = delete;

private:

};

class MmffAtoms {
public:

	static const MmffAtoms & getInstance();

	virtual ~MmffAtoms() {
	}

	MmffAtoms(const MmffAtoms & rhs) = delete;
	MmffAtoms & operator =(const MmffAtoms & rhs) = delete;
	MmffAtoms(MmffAtoms && rhs) = delete;
	MmffAtoms & operator =(MmffAtoms && rhs) = delete;

private:

	MmffAtoms() {
	}

	static std::vector<unique_ptr<MmffAtom>> & mmffAtomDefinitions();
};

class MmffAtomSmarts {
public:
	MmffAtomSmarts(const string & sym, const int no, const int mmffNo,
			const int h, const int charge, const int v, const string & n, const string & smi) :
			atomSymbol(sym), atomicNo(no), mmff94No(mmffNo), hType(h), formalCharge12(
					charge), valance(v), name(n), smiles(smi) {
	}

	virtual ~MmffAtomSmarts() {
	}

	MmffAtomSmarts(const MmffAtomSmarts & rhs) = delete;
	MmffAtomSmarts & operator =(const MmffAtomSmarts & rhs) = delete;
	MmffAtomSmarts(MmffAtomSmarts && rhs) = delete;
	MmffAtomSmarts & operator =(MmffAtomSmarts && rhs) = delete;

    const int getAtomicNo() const {
        return atomicNo;
    }

    const string& getAtomSymbol() const {
        return atomSymbol;
    }

    const int getFormalCharge12() const {
        return formalCharge12;
    }

    const int getType() const {
        return hType;
    }

    const int getMmff94No() const {
        return mmff94No;
    }

    const string& getName() const {
        return name;
    }

    const string& getSmiles() const {
        return smiles;
    }

    const int getValance() const {
        return valance;
    }

private:
	const string atomSymbol;
	const int atomicNo, mmff94No, hType;
	const int formalCharge12, valance;
	const string name, smiles;
};

} /* namespace GarethFF */

#endif /* SRC_FF_MMFFATOMS_H_ */
