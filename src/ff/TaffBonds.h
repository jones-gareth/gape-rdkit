/*
 * TaffBonds.h
 *
 *  Created on: Mar 12, 2016
 *      Author: gjones
 */

#ifndef SRC_FF_TAFFBONDS_H_
#define SRC_FF_TAFFBONDS_H_

#include <map>
#include <memory>
#include <vector>

#include "../mol/AtomType.h"
#include "../mol/BondType.h"
#include "../mol/Molecule.h"

namespace GarethMol {
class Molecule;
} /* namespace GarethMol */

namespace GarethFF {

using namespace GarethMol;
using namespace std;

class TaffBond {
    using AtomTypeId = GarethMol::AtomType::AtomTypeId;
    using BondTypeId = GarethMol::BondType::BondTypeId;

public:

    TaffBond(const AtomTypeId a1, const AtomTypeId a2, const BondTypeId b,
            const double l, const double k_) :
            atom1TypeId(a1), atom2TypeId(a2), bondTypeId(b), length(l), k(k_) {
        weight = 2;
        if (atom1TypeId == AtomTypeId::WILD ) {
            weight--;
        }
        if (atom2TypeId == AtomTypeId::WILD) {
            weight--;
        }
    }

    virtual ~TaffBond() {
    }

    TaffBond(const TaffBond & rhs) = delete;
    TaffBond & operator =(const TaffBond & rhs) = delete;
    TaffBond(TaffBond && rhs) = delete;
    TaffBond & operator =(TaffBond && rhs) = delete;

    const AtomTypeId getAtom1TypeId() const {
        return atom1TypeId;
    }

    const AtomTypeId getAtom2TypeId() const {
        return atom2TypeId;
    }

    const BondTypeId getBondTypeId() const {
        return bondTypeId;
    }

    const double getK() const {
        return k;
    }

    const double getLength() const {
        return length;
    }

    int getWeight() const {
        return weight;
    }

private:
    const AtomTypeId atom1TypeId, atom2TypeId;
    const BondTypeId bondTypeId;
    const double length;
    const double k;
    int weight;

};

class TaffBonds {
public:

    static const TaffBonds & getInstance();

    virtual ~TaffBonds() {
    }

    TaffBonds(const TaffBonds & rhs) = delete;
    TaffBonds & operator =(const TaffBonds & rhs) = delete;
    TaffBonds(TaffBonds && rhs) = delete;
    TaffBonds & operator =(TaffBonds && rhs) = delete;

    map<int, TaffBond *> parametersForMolecule(const Molecule &mol) const;

private:

    TaffBonds() {
    }

    static std::vector<unique_ptr<TaffBond>> & taffBondDefinitions();

    static  unique_ptr<TaffBond>  & defaultParameters();
};
} /* namespace GarethFF */

#endif /* SRC_FF_TAFFBONDS_H_ */
