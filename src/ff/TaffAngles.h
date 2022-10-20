/*
 * TaffAngle.h
 *
 *  Created on: Mar 12, 2016
 *      Author: gjones
 */

#ifndef SRC_FF_TAFFANGLES_H_
#define SRC_FF_TAFFANGLES_H_

#include <memory>
#include "../mol/Atom.h"

namespace GarethFF {

using namespace GarethMol;
using namespace std;

/**
 * A class to represent angle parameters for the Tripos force field
 */
class TaffAngle {

public:
	  using AtomTypeId = GarethMol::AtomType::AtomTypeId;

	  TaffAngle(const AtomTypeId a1, const AtomTypeId a2, const AtomTypeId a3,
	            const double a, const double k_) :
	            atom1TypeId(a1), atom2TypeId(a2), atom3TypeId(a3), angle(a), k(k_) {
	        weight = 3;
	        if (atom1TypeId == AtomTypeId::WILD) {
	            weight--;
	        }
	        if (atom2TypeId == AtomTypeId::WILD) {
	            weight--;
	        }
	        if (atom3TypeId == AtomTypeId::WILD) {
	            weight--;
	        }
	    }

	  virtual ~TaffAngle() {
    }

    TaffAngle(const TaffAngle & rhs) = delete;
    TaffAngle & operator =(const TaffAngle & rhs) = delete;
    TaffAngle(TaffAngle && rhs) = delete;
    TaffAngle & operator =(TaffAngle && rhs) = delete;

    const double getAngle() const {
        return angle;
    }

    const AtomTypeId getAtom1TypeId() const {
        return atom1TypeId;
    }

    const AtomTypeId getAtom2TypeId() const {
        return atom2TypeId;
    }

    const AtomTypeId getAtom3TypeId() const {
        return atom3TypeId;
    }

    const double getK() const {
        return k;
    }

    const int getWeight() const {
        return weight;
    }

private:
    const AtomTypeId atom1TypeId, atom2TypeId, atom3TypeId;
    const double angle;
    const double k;
    int weight;


};

/**
 * A singleton class to hold all angle parameters for the Tripos forcefield and
 * handle parameterization of a molecule.
 */
class TaffAngles {
public:

    static TaffAngles & getInstance();

    virtual ~TaffAngles() {
    }

    TaffAngles(const TaffAngles & rhs) = delete;
    TaffAngles & operator =(const TaffAngles & rhs) = delete;
    TaffAngles(TaffAngles && rhs) = delete;
    TaffAngles & operator =(TaffAngles && rhs) = delete;

    map<int, TaffAngle *> parametersForMolecule(const Molecule & mol) const;

private:
    TaffAngles() {
    }

    static  std::vector<unique_ptr<TaffAngle>> & taffAngleDefinitions();
    static  map<AtomType::Geometry, unique_ptr<TaffAngle>>  & defaultParameters();

};
} /* namespace GarethFF */

#endif /* SRC_FF_TAFFANGLES_H_ */
