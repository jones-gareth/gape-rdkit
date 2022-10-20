/*
 * TaffTorsions.h
 *
 *  Created on: Mar 15, 2016
 *      Author: gareth
 */

#ifndef SRC_FF_TAFFTORSIONS_H_
#define SRC_FF_TAFFTORSIONS_H_

#include <memory>
#include "../mol/Atom.h"
#include "../mol/Molecule.h"

namespace GarethFF {

using namespace GarethMol;
using namespace std;

class TaffTorsion {

    using AtomTypeId = GarethMol::AtomType::AtomTypeId;
    using BondTypeId = GarethMol::BondType::BondTypeId;

public:
    TaffTorsion(const AtomTypeId a1, const AtomTypeId a2, const AtomTypeId a3,
            const AtomTypeId a4, const BondTypeId bt, const double k_,
            const double p_) :
            atom1TypeId(a1), atom2TypeId(a2), atom3TypeId(a3), atom4TypeId(a4), bondTypeId(
                    bt), k(k_), p(p_) {
        weight = 4;
        if (atom1TypeId == AtomTypeId::WILD) {
            weight--;
        }
        if (atom2TypeId == AtomTypeId::WILD) {
            weight--;
        }
        if (atom3TypeId == AtomTypeId::WILD) {
            weight--;
        }
        if (atom4TypeId == AtomTypeId::WILD) {
            weight--;
        }
    }

    virtual ~TaffTorsion() {
    }

    TaffTorsion(const TaffTorsion & rhs) = delete;
    TaffTorsion & operator =(const TaffTorsion & rhs) = delete;
    TaffTorsion(TaffTorsion && rhs) = delete;
    TaffTorsion & operator =(TaffTorsion && rhs) = delete;

    const AtomTypeId getAtom1TypeId() const {
        return atom1TypeId;
    }

    const AtomTypeId getAtom2TypeId() const {
        return atom2TypeId;
    }

    const AtomTypeId getAtom3TypeId() const {
        return atom3TypeId;
    }

    const AtomTypeId getAtom4TypeId() const {
        return atom4TypeId;
    }

    const double getK() const {
        return k;
    }

    const double getP() const {
        return p;
    }

    int getWeight() const {
        return weight;
    }

    const BondTypeId getBondTypeId() const {
        return bondTypeId;
    }

private:
    const AtomTypeId atom1TypeId, atom2TypeId, atom3TypeId, atom4TypeId;
    const BondTypeId bondTypeId;
    const double k;
    const double p;
    int weight;
};

class MatchedTaffTorsion {
public:
    MatchedTaffTorsion(bool rev, const TaffTorsion * const taffTor,
            const Torsion * const tor) :
            reverse(rev), taffTorsion(taffTor), torsion(tor) {
    }

    virtual ~MatchedTaffTorsion() {
    }

    MatchedTaffTorsion(const MatchedTaffTorsion & rhs) = delete;
    MatchedTaffTorsion & operator =(const MatchedTaffTorsion & rhs) = delete;
    MatchedTaffTorsion(MatchedTaffTorsion && rhs) = delete;
    MatchedTaffTorsion & operator =(MatchedTaffTorsion && rhs) = delete;

    const bool isReverse() const {
        return reverse;
    }

    const TaffTorsion* const & getTaffTorsion() const {
        return taffTorsion;
    }

    const Torsion* const & getTorsion() const {
        return torsion;
    }

private:
    const bool reverse;
    const TaffTorsion * const taffTorsion;
    const Torsion * const torsion;
};

class TaffTorsions {
public:

    static const TaffTorsions & getInstance();

    virtual ~TaffTorsions() {
    }

    TaffTorsions(const TaffTorsions & rhs) = delete;
    TaffTorsions & operator =(const TaffTorsions & rhs) = delete;
    TaffTorsions(TaffTorsions && rhs) = delete;
    TaffTorsions & operator =(TaffTorsions && rhs) = delete;

    map<int, unique_ptr<MatchedTaffTorsion>> parametersForMolecule(
            const Molecule & mol) const;
private:
    TaffTorsions() {
    }

    static std::vector<unique_ptr<TaffTorsion>> & taffTorsionDefinitions();

};

} /* namespace GarethFF */

#endif /* SRC_FF_TAFFTORSIONS_H_ */
