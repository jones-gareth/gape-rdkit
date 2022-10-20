/*
 * Taff.h
 *
 *  Created on: Mar 8, 2016
 *      Author: gareth
 */

#ifndef SRC_FF_TAFF_H_
#define SRC_FF_TAFF_H_

#include "../mol/Molecule.h"
#include "../util/Array2D.h"
#include "TaffAtoms.h"
#include "TaffAngles.h"
#include "TaffBonds.h"
#include "TaffOops.h"
#include "TaffTorsions.h"

namespace GarethFF {

using namespace GarethMol;
using namespace std;

using AtomTypeId = AtomType::AtomTypeId;
using BondTypeId = BondType::BondTypeId;

/**
 * A class to implement the Tripos Associates Forcefield
 */
class Taff {
public:
    /**
     * Constructor- initialize using a molecule
     *
     * @param mol
     */
    Taff(const Molecule & mol) :
            molecule(mol) {
        setupMolecule();
    }

    virtual ~Taff() {
    }

    Taff(const Taff & rhs) = delete;
    Taff & operator =(const Taff & rhs) = delete;
    Taff(Taff && rhs) = delete;
    Taff & operator =(Taff && rhs) = delete;

    /**
     * Determine the total VDW energy for the molecule with the given coordinates
     * @param coords
     * @return
     */
    double vdwEnergy(const CoordMatrix & coords) const;

    /**
     * Determine the total angle bend energy for the molecule with the given coordinates
     *
     * @param coords
     * @return
     */
    double angleEnergy(const CoordMatrix & coords) const;

    /**
     * Determine the total bond stretch energy for the molecule with the given coordinates
     *
     * @param coords
     * @return
     */
    double bondEnergy(const CoordMatrix & coords) const;

    /**
     * Determine the total out of plane energy for the molecule with the given co-ordinates
     * @param coords
     * @return
     */
    double oopEnergy(const CoordMatrix & coords) const;

    /**
     * Determine the total torsional energy for the molecule with the given co-ordinates
     * @param coords
     * @return
     */
    double torsionEnergy(const CoordMatrix & coords) const;

    /**
     * Determines the total energy for the molecule with the coordinates
     *
     * @param coords
     * @return
     */
    double energy(const CoordMatrix & coords);

    /**
     *
     * @param atom1No
     * @param atom2No
     * @param coords
     * @return
     */
    double vdwPairEnergy(int atom1No, int atom2No,
            const CoordMatrix & coords) const;

    /**
     * Determines the angle energy for the given angle for the molecule using the
     * supplied coordinates.
     *
     * @param angleNo
     * @param coords
     * @return
     */
    double angleEnergy(int angleNo, const CoordMatrix & coords) const;

    /**
     * Determines the torsion energy for the given angle for the molecule using the
     * supplied coordinates.
     *
     * @param torsionNo
     * @param coords
     * @return
     */
    double torsionEnergy(int torsionNo, const CoordMatrix & coords) const;

    /**
     * Determines the bond energy for the given angle for the molecule using the
     * supplied coordinates.
     *
     * @param bondNo
     * @param coords
     * @return
     */
    double bondEnergy(int bondNo, const CoordMatrix & coords) const;

    /**
     * Determines the out of plane energy for the given atom for the molecule using the
     * supplied coordinates.
     *
     * @param atomNo
     * @param coords
     * @return
     */
    double oopEnergy(int atomNo, const CoordMatrix & coords) const;

    /**
     * Matches a forcefield atom type to a molecule atom type
     * @param moleculeAtomType
     * @param ffAtomType
     * @return
     */
    static bool matchAtom(AtomTypeId moleculeAtomType, AtomTypeId ffAtomType);

    /**
     * Matches a molecule bond type to a forcefield bond type
     * @param molecule
     * @param moleculeBond
     * @param ffBondType
     * @return
     */
    static bool matchBond(const Molecule & molecule, const Bond & moleculeBond,
            const BondTypeId ffBondType);

    /**
     * Determines the sum of all contributions to the force field that an atom
     * makes using the supplied molecular coordinates
     *
     * @param atom
     * @param coords
     * @return
     */
    double atomEnergy(const Atom & atom, const CoordMatrix & coords);

    /**
     * Determines the sum of all VDW interactions for this atom using the supplied
     * molecular coordinates
     *
     * @param atom
     * @param coords
     * @return
     */
    double atomVdwEnergy(const Atom & atom, const CoordMatrix & coords) const;

    /**
     * Determines the sum of all bond stretch interactions for this atom using the supplied
     * molecular coordinates
     *
     * @param atom
     * @param coords
     * @return
     */
    double atomBondEnergy(const Atom & atom, const CoordMatrix & coords) const;

    /**
     * Determines the sum of all bond angle pairs that contain this atom using the supplied
     * molecular coordinates
     *
     * @param atom
     * @param coords
     * @return
     */
    double atomAngleEnergy(const Atom & atom, const CoordMatrix & coords) const;

    /**
     * Determines the sum of all torsion energies that contain this atom using the supplied
     * molecular coordinates
     *
     * @param atom
     * @param coords
     * @return
     */
    double atomTorsionEnergy(const Atom & atom,
            const CoordMatrix & coords) const;

    /**
     * Minimize a molecule by performing simplex minimization on each atom
     * @param coords
     * @param nCycles
     * @return
     */
    double minimizeAtoms(CoordMatrix & coords, const int nCycles = 100,
            const double vdwDistanceCutoff = .0, const double step = 0.25);

    /**
     * Minimize a molecule by performing simplex minimization on all atoms
     * simultaneously
     *
     * @param coords
     * @param nCycles
     * @return
     */
    double minimizeMolecule(CoordMatrix & coords, const int nCycles = 5000,
            const double step = 0.25);

    // getters

    const Molecule &getMolecule() const {
        return molecule;
    }

    double getEVdw() const {
        return eVdw;
    }

    double getEBond() const {
        return eBond;
    }

    double getEAngle() const {
        return eAngle;
    }

    double getETorsion() const {
        return eTorsion;
    }

    double getEOop() const {
        return eOop;
    }

    double getAtomEVdw() const {
        return atomEVdw;
    }

    double getAtomEBond() const {
        return atomEBond;
    }

    double getAtomEAngle() const {
        return atomEAngle;
    }

    double getAtomETorsion() const {
        return atomETorsion;
    }

    double getAtomEOop() const {
        return atomEOop;
    }

    void setVdwMaxDistance(double vdwMaxDistance) {
        this->vdwMaxSqrDistance = vdwMaxDistance * vdwMaxDistance;
    }

private:
    const Molecule & molecule;
    map<int, TaffAtom *> taffAtoms;
    map<int, TaffAngle *> taffAngles;
    map<int, TaffBond *> taffBonds;
    map<int, TaffOop *> taffOops;
    map<int, unique_ptr<MatchedTaffTorsion>> taffTorsions;
    unique_ptr<Array2D<bool>> pairsToCheck;
    double vdwMaxSqrDistance = .0;

    void setupMolecule();

    double eVdw, eBond, eAngle, eTorsion, eOop;

    double atomEVdw, atomEBond, atomEAngle, atomETorsion, atomEOop;

};

} /* namespace GarethFF */

#endif /* SRC_FF_TAFF_H_ */
