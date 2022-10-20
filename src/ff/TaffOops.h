/*
 * TaffOops.h
 *
 *  Created on: Mar 13, 2016
 *      Author: gjones
 */

#ifndef SRC_FF_TAFFOOPS_H_
#define SRC_FF_TAFFOOPS_H_

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

class TaffOop {
    using AtomTypeId = GarethMol::AtomType::AtomTypeId;

public:
    TaffOop( const double k_, const AtomTypeId id_) :
         k(k_), atomTypeId(id_) {
    }

    virtual ~TaffOop() {
    }

    TaffOop(const TaffOop & rhs) = delete;
    TaffOop & operator =(const TaffOop & rhs) = delete;
    TaffOop(TaffOop && rhs) = delete;
    TaffOop & operator =(TaffOop && rhs) = delete;

    const AtomTypeId getAtomTypeId() const {
        return atomTypeId;
    }

    const double getK() const {
        return k;
    }

private:
    const double k;
    const AtomTypeId atomTypeId;

};

class TaffOops {
public:

    static const TaffOops & getInstance();

    virtual ~TaffOops() {
    }

    TaffOops(const TaffOops & rhs) = delete;
    TaffOops & operator =(const TaffOops & rhs) = delete;
    TaffOops(TaffOops && rhs) = delete;
    TaffOops & operator =(TaffOops && rhs) = delete;

    map<int, TaffOop *> parametersForMolecule(const Molecule & mol) const;

private:

    TaffOops() {
    }

    static std::vector<unique_ptr<TaffOop>> & taffOopDefinitions() ;
};
} /* namespace GarethFF */

#endif /* SRC_FF_TAFFOOPS_H_ */
