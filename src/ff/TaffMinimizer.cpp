/*
 * TaffMinimizer.cpp
 *
 *  Created on: May 6, 2016
 *      Author: gjones
 */

#include "TaffMinimizer.h"
#include "../util/Reporter.h"

using namespace GarethUtil;
using namespace Eigen;

namespace GarethFF {

double TaffSimplexAtom::minimize(Atom & atom, CoordMatrix & coords) {
    this->coords = &coords;
    this->atom = &atom;

    Simplex<TaffSimplexAtom> simplex(*this);
    double energy = simplex.minimize();
    copyToCoordinates(values);
    return energy;
}

int TaffSimplexAtom::getProblemSize() {
    return 3;
}

void TaffSimplexAtom::copyToCoordinates(const double vals[]) {
    auto atomNo = atom->getAtomNo();
    auto atomCoords = coords->col(atomNo);
    for (auto i = 0; i < 3; i++) {
        atomCoords[i] = vals[i];
    }
}

double * TaffSimplexAtom::getStart() {
    auto atomNo = atom->getAtomNo();
    auto atomCoords = coords->col(atomNo);
    for (auto i = 0; i < 3; i++) {
        values[i] = atomCoords[i];
    }
    return values;
}

double TaffSimplexAtom::evaluate(double values[]) {
    copyToCoordinates(values);
    double energy = taff.atomEnergy(*atom, *coords);
    return energy;
}

double * TaffSimplexAtom::getSteps() {
    for (auto i=0; i<3; i++)
        steps[i] = step;
    return steps;
}

double TaffSimplexAtom::getConvergenceCriteria() {
    return 1e-2;
}

int TaffSimplexAtom::checkConvergenceIterations() {
    return 10;
}

int TaffSimplexAtom::maxIterations() {
    return 100;
}

double TaffSimplexMolecule::minimize(CoordMatrix & coords) {
    this->coords = &coords;

    Simplex<TaffSimplexMolecule> simplex(*this);
    double energy = simplex.minimize();
    coords.block(0, 0, 3, molecule.nAtoms()) = atomCoords;
    REPORT(Reporter::DEBUG) << "Final simplex molecule energy " << energy;
    return energy;
}

int TaffSimplexMolecule::getProblemSize() {
    auto nAtoms = molecule.nAtoms();
    return 3*nAtoms;
}

double * TaffSimplexMolecule::getStart() {
    auto nAtoms = molecule.nAtoms();
    atomCoords = coords->block(0, 0, 3, nAtoms);
    return atomCoords.data();
}

double TaffSimplexMolecule::evaluate(double values[]) {
    auto nAtoms = molecule.nAtoms();
    Map<const MatrixXd> mf(values, 3, nAtoms);
    coords->block(0, 0, 3, nAtoms) = mf;
    double energy = taff.energy(*coords);
    REPORT(Reporter::DETAIL) << "Simplex molecule energy evaluation " << energy;
    return energy;
}

double * TaffSimplexMolecule::getSteps() {
    auto size = molecule.nAtoms() * 3;
    steps = new double[size];
    for (auto i = 0ul; i < size; i++) {
        steps[i] = step;
    }
    return steps;
}

double TaffSimplexMolecule::getConvergenceCriteria() {
    return 1e-3;
}

int TaffSimplexMolecule::checkConvergenceIterations() {
    return 10;
}

int TaffSimplexMolecule::maxIterations() {
    return nCycles;
}

} /* namespace GarethFF */
