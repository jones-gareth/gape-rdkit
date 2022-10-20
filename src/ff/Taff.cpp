/*
 rf  * Taff.cpp
 *
 *  Created on: Mar 8, 2016
 *      Author: gareth
 */

#include "Taff.h"

#include <boost/range/irange.hpp>
#include <Eigen/Dense>
#include <cassert>
#include <cmath>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "../mol/Atom.h"
#include "../mol/AtomType.h"
#include "../mol/Bond.h"
#include "../util/CoordOps.h"
#include "../util/Reporter.h"
#include "TaffMinimizer.h"

namespace GarethFF {

/**
 *
 * Use this matching function for bond definitions instead of the atom type one.
 * I need to investigate the Tripos TAFF_BOND_STRETCH file more- it seems like you
 * don't do matching for pairs like "HEV H", as you get crappy bonds and the more
 *  reasonable "* H" is never matched.  Need to chek out the reference and see if
 *  it is any clearer
 *
 * @param moleculeAtomType
 * @param ffAtomType
 * @return
 */
bool Taff::matchAtom(const AtomTypeId moleculeAtomType,
        const AtomTypeId ffAtomType) {
    assert(moleculeAtomType != AtomTypeId::WILD);
    assert(ffAtomType != AtomTypeId::OAR);

    if (ffAtomType == AtomTypeId::WILD) {
        return true;
    }
    // O.ar is not in any of the TAFF forcefield files- use O.3 as a surrogate
    else if (moleculeAtomType == AtomTypeId::OAR) {
        return AtomTypeId::O3 == ffAtomType;
    } else {
        return moleculeAtomType == ffAtomType;
    }
}

bool Taff::matchBond(const Molecule & molecule, const Bond & moleculeBond,
        const BondTypeId ffBondType) {
    return moleculeBond.sybylBondType(molecule) == ffBondType;
}

void Taff::setupMolecule() {
    taffAtoms = TaffAtoms::getInstance().parametersForMolecule(molecule);
    taffAngles = TaffAngles::getInstance().parametersForMolecule(molecule);
    taffBonds = TaffBonds::getInstance().parametersForMolecule(molecule);
    taffOops = TaffOops::getInstance().parametersForMolecule(molecule);
    taffTorsions = TaffTorsions::getInstance().parametersForMolecule(molecule);

    int nAtoms = molecule.nAtoms();
    pairsToCheck = make_unique<Array2D<bool>>(molecule.nAtoms(), nAtoms);
    for (auto i = 0; i < nAtoms; i++) {
        for (auto j = i + 1; j < nAtoms; j++) {
            auto check = true;
            if (!AtomType::isRealAtom(molecule.getAtom(i).getAtomTypeId()))
                check = false;
            if (!AtomType::isRealAtom(molecule.getAtom(j).getAtomTypeId()))
                check = false;
            if (molecule.is12Bonded(i, j)) {
                check = false;
            }
            if (molecule.is13Bonded(i, j)) {
                check = false;
            }
            pairsToCheck->set(i, j, check);
            pairsToCheck->set(j, i, check);

            if (check) {
                if (taffAtoms.count(i) == 0) {
                    REPORT(Reporter::WARN)
                            << "Missing VDW parameters for atom "
                            << to_string(i + 1);
                    assert(false);
                }
                if (taffAtoms.count(j) == 0) {
                    REPORT(Reporter::WARN)
                            << "Missing VDW parameters for atom "
                            << to_string(j + 1);
                    assert(false);
                }
            }
        }
    }
}

double Taff::vdwEnergy(const CoordMatrix & coords) const {
    auto nAtoms = molecule.nAtoms();
    double eVdw = .0;

    for (auto i = 0ul; i < nAtoms; i++) {
        for (auto j = i + 1; j < nAtoms; j++) {
            if (pairsToCheck->get(i, j))
                eVdw += vdwPairEnergy(i, j, coords);
        }
    }

    REPORT(Reporter::TRACE) << "Molecular VDW energy " << eVdw;
    return eVdw;
}

double Taff::vdwPairEnergy(int i, int j, const CoordMatrix & coords) const {
    if (!pairsToCheck->get(i, j)) {
        return .0;
    }

    auto taffAtom1 = taffAtoms.at(i);
    auto taffAtom2 = taffAtoms.at(j);
    auto sqrD = sqrDistance(coords.col(i), coords.col(j));
    if (vdwMaxSqrDistance > .0 && sqrD > vdwMaxSqrDistance) {
        return .0;
    }

    auto r = taffAtom1->getR() + taffAtom2->getR();
    sqrD = sqrD / (r * r);
    auto pow6 = sqrD * sqrD * sqrD;
    auto pow12 = pow6 * pow6;
    auto kMean = sqrt(taffAtom1->getK() * taffAtom2->getK());
    auto atomVdw = kMean * (1.0 / pow12 - 2.0 / pow6);

    REPORT(Reporter::TRACE) << "VDW energy between atoms " << to_string(i + 1)
            << " and " << to_string(j + 1) << ": " << atomVdw;
    return atomVdw;
}

double Taff::angleEnergy(const CoordMatrix & coords) const {
    double eAngle = 0;
    for (auto angleNo : boost::irange(0ul, molecule.getAngles().size())) {
        eAngle += angleEnergy(angleNo, coords);
    }

    REPORT(Reporter::TRACE) << "Molecular Angle energy " << eAngle;
    return eAngle;
}

double Taff::angleEnergy(int angleNo, const CoordMatrix & coords) const {
    const auto & molAngle = molecule.getAngles().at(angleNo);
    double angleRadians = angle(coords.col(molAngle->getAtom1()),
            coords.col(molAngle->getAtom2()), coords.col(molAngle->getAtom3()));
    double ang = angleRadians * 180.0 / M_PI;

    const auto & taffAngle = taffAngles.at(angleNo);
    auto diff = taffAngle->getAngle() - ang;
    REPORT(Reporter::TRACE) << "Angle is " << ang << " TAFF angle "
            << taffAngle->getAngle() << " k " << taffAngle->getK();
    auto e = taffAngle->getK() * diff * diff;
    REPORT(Reporter::TRACE) << "Angle energy between atoms "
            << to_string(molAngle->getAtom1() + 1) << ", "
            << to_string(molAngle->getAtom2() + 1) << " and "
            << to_string(molAngle->getAtom3() + 1) << ": " << e;
    return e;
}

double Taff::bondEnergy(const CoordMatrix & coords) const {
    double eBond = 0;
    for (auto bondNo : boost::irange(0ul, molecule.nBonds())) {
        eBond += bondEnergy(bondNo, coords);
    }

    REPORT(Reporter::TRACE) << "Molecular Bond energy " << eBond;
    return eBond;
}

double Taff::bondEnergy(int bondNo, const CoordMatrix & coords) const {
    const auto & bond = molecule.getBond(bondNo);
    if (!AtomType::isRealAtom(bond.getAtom1().getAtomTypeId())) {
        return .0;
    }
    if (!AtomType::isRealAtom(bond.getAtom2().getAtomTypeId())) {
        return .0;
    }

    auto atom1No = bond.getAtom1().getAtomNo();
    auto atom2No = bond.getAtom2().getAtomNo();
    const auto & taffBond = taffBonds.at(bondNo);
    auto len = distance(coords.col(atom1No), coords.col(atom2No));
    auto d = len - taffBond->getLength();
    REPORT(Reporter::TRACE) << "Bond length is " << len << " Taff length "
            << taffBond->getLength() << " K " << taffBond->getK();
    auto e = taffBond->getK() * (d * d);
    REPORT(Reporter::TRACE) << "Bond energy between atoms "
            << to_string(atom1No + 1) << " and " << to_string(atom2No + 1)
            << " is " << e;
    return e;
}

double Taff::oopEnergy(const CoordMatrix & coords) const {
    auto eOop = .0;
    for (auto i : boost::irange(0ul, molecule.nAtoms())) {
        eOop += oopEnergy(i, coords);
    }

    REPORT(Reporter::TRACE) << "Molecular Out of plane  energy " << eOop;
    return eOop;
}

double Taff::oopEnergy(int atomNo, const CoordMatrix & coords) const {
    if (taffOops.count(atomNo) == 0) {
        return .0;
    }
    auto neighbours =
            molecule.getAtom(atomNo).getNeighbourhood().getAtomNeighbours();
    assert(neighbours.size() == 3);
    auto height = outOfPlaneHeight(coords.col(atomNo),
            coords.col(neighbours.at(0)->getAtomNo()),
            coords.col(neighbours.at(1)->getAtomNo()),
            coords.col(neighbours.at(2)->getAtomNo()));

    auto k = taffOops.at(atomNo)->getK();
    auto atomOop = k * height * height;
    REPORT(Reporter::TRACE) << "K " << k << " height " << height;
    REPORT(Reporter::TRACE) << "Atom Out of plane energy for "
            << to_string(atomNo + 1) << " : " << atomOop;
    return atomOop;
}

double Taff::torsionEnergy(const CoordMatrix & coords) const {
    auto eTorsion = .0;
    for (auto torsionNo : boost::irange(0ul, molecule.getTorsions().size())) {
        eTorsion += torsionEnergy(torsionNo, coords);
    }
    REPORT(Reporter::TRACE) << "Molecular torsion  energy " << eTorsion;
    return eTorsion;
}

double Taff::torsionEnergy(int torsionNo, const CoordMatrix & coords) const {
    if (taffTorsions.count(torsionNo) == 0) {
        return .0;
    }
    const auto & matchedTorsion = taffTorsions.at(torsionNo);
    const auto & torsion = molecule.getTorsions().at(torsionNo);
    assert(matchedTorsion->getTorsion() == torsion.get());

    auto a = coords.col(torsion->getAtom1());
    auto b = coords.col(torsion->getAtom2());
    auto c = coords.col(torsion->getAtom3());
    auto d = coords.col(torsion->getAtom4());

    double angle =
            matchedTorsion->isReverse() ?
                    torsionAngle(d, c, b, a) : torsionAngle(a, b, c, d);

    auto p = matchedTorsion->getTaffTorsion()->getP();
    auto k = matchedTorsion->getTaffTorsion()->getK();
    auto sign = (p > 0) ? 1.0 : -1.0;
    auto e = k * (1 + sign * cos(p * angle));

    REPORT(Reporter::TRACE) << "Angle is " << angle << " TAFF torsion P " << p
            << " k " << k;
    REPORT(Reporter::TRACE) << "Torsion energy between atoms "
            << to_string(torsion->getAtom1() + 1) << ", "
            << to_string(torsion->getAtom2() + 1) << ", "
            << to_string(torsion->getAtom3() + 1) << " and "
            << to_string(torsion->getAtom4() + 1) << ": " << e;

    return e;
}

double Taff::energy(const CoordMatrix & coords) {
    eVdw = vdwEnergy(coords);
    eAngle = angleEnergy(coords);
    eBond = bondEnergy(coords);
    eTorsion = torsionEnergy(coords);
    eOop = oopEnergy(coords);
    auto energy = eVdw + eAngle + eBond + eTorsion + eOop;

    REPORT(Reporter::DEBUG) << "Total energy " << energy << " eVdw " << eVdw
            << " eAngle " << eAngle << " eBond " << eBond << " eTorsion "
            << eTorsion << " eOop " << eOop;

    return energy;
}

double Taff::atomEnergy(const Atom & atom, const CoordMatrix & coords) {
    atomEVdw = atomVdwEnergy(atom, coords);
    atomEAngle = atomAngleEnergy(atom, coords);
    atomEBond = atomBondEnergy(atom, coords);
    atomETorsion = atomTorsionEnergy(atom, coords);
    atomEOop = oopEnergy(atom.getAtomNo(), coords);

    auto energy = atomEVdw + atomEAngle + atomEBond + atomETorsion + atomEOop;

    REPORT(Reporter::TRACE) << "Atom energy for atom "
            << to_string(atom.getAtomNo() + 1) << ": " << energy << " eVdw "
            << atomEVdw << " eAngle " << atomEAngle << " eBond " << atomEBond
            << " eTorsion " << atomETorsion << " eOop " << atomEOop;

    return energy;
}

double Taff::atomVdwEnergy(const Atom & atom,
        const CoordMatrix & coords) const {

    auto nAtoms = molecule.nAtoms();
    auto atomNo = atom.getAtomNo();
    double eVdw = .0;

    for (auto i = 0ul; i < nAtoms; i++) {
        eVdw += vdwPairEnergy(i, atomNo, coords);
    }

    REPORT(Reporter::TRACE) << "Atom VDW energy " << eVdw;
    return eVdw;
}

double Taff::atomBondEnergy(const Atom & atom,
        const CoordMatrix & coords) const {
    double eBond = .0;
    for (const auto bond : atom.getNeighbourhood().getBonds()) {
        eBond += bondEnergy(bond->getBondNo(), coords);
    }
    REPORT(Reporter::TRACE) << "Atom bond energy " << eBond;
    return eBond;
}

double Taff::atomAngleEnergy(const Atom & atom,
        const CoordMatrix & coords) const {
    double eAngle = .0;
    for (const auto angle : atom.getNeighbourhood().getAngles()) {
        eAngle += angleEnergy(angle->getAngleNo(), coords);
    }
    REPORT(Reporter::TRACE) << "Atom angle energy " << eBond;
    return eAngle;
}

double Taff::atomTorsionEnergy(const Atom & atom,
        const CoordMatrix & coords) const {
    double eTorsion = .0;
    for (const auto torsion : atom.getNeighbourhood().getTorsions()) {
        eTorsion += torsionEnergy(torsion->getTorsionNo(), coords);
    }
    REPORT(Reporter::TRACE) << "Atom torsion energy " << eBond;
    return eTorsion;
}

double Taff::minimizeAtoms(CoordMatrix & coords, const int nCycles,
        const double vdwDistanceCutoff, const double step) {
    auto molEnergy = energy(coords);
    REPORT(Reporter::DETAIL) << "Starting energy " << molEnergy;
    TaffSimplexAtom atomMinimizer(*this);
    atomMinimizer.setStep(step);
    for (auto no = 0; no < nCycles; no++) {
        setVdwMaxDistance(vdwDistanceCutoff);
        for (auto & atom : molecule.getAtoms()) {
            atomMinimizer.minimize(*atom, coords);
        }
        setVdwMaxDistance(0);
        molEnergy = energy(coords);
        REPORT(Reporter::DETAIL) << "Cycle " << to_string(no + 1) << " energy "
                << molEnergy;
    }
    REPORT(Reporter::DETAIL) << "Ending energy " << molEnergy;
    return molEnergy;
}

double Taff::minimizeMolecule(CoordMatrix & coords, const int nCycles,
        const double step) {
    auto molEnergy = energy(coords);
    REPORT(Reporter::DETAIL) << "Starting energy " << molEnergy;
    TaffSimplexMolecule molMinimizer(*this, nCycles);
    molMinimizer.setStep(step);
    molEnergy = molMinimizer.minimize(coords);
    REPORT(Reporter::DETAIL) << "Ending energy " << molEnergy;
    return molEnergy;
}

} /* namespace GarethFF */
