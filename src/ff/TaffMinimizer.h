/*
 * TaffMinimizer.h
 *
 *  Created on: May 6, 2016
 *      Author: gjones
 */

#ifndef SRC_FF_TAFFMINIMIZER_H_
#define SRC_FF_TAFFMINIMIZER_H_

#include "Taff.h"
#include "../util/Simplex.h"

using namespace GarethUtil;

namespace GarethFF {

/**
 * A class to handle simplex minimization of a single atom
 *
 */
class TaffSimplexAtom {
    friend class Simplex<TaffSimplexAtom> ;

public:
    TaffSimplexAtom(Taff & taff_) :
            taff(taff_) {
    }

    virtual ~TaffSimplexAtom() {
    }

    TaffSimplexAtom(const TaffSimplexAtom & rhs) = delete;
    TaffSimplexAtom & operator =(const TaffSimplexAtom & rhs) = delete;
    TaffSimplexAtom(TaffSimplexAtom && rhs) = delete;
    TaffSimplexAtom & operator =(TaffSimplexAtom && rhs) = delete;

    /**
     * Minimize the atom using the supplied molecular coordinates
     *
     * @param atom
     * @param coords
     * @return
     */
    double minimize(Atom & atom, CoordMatrix & coords);

    void setStep(double step) {
          this->step = step;
      }

private:
    /**
     * Forcefield parameters for molecule
     */
    Taff & taff;

    /**
     * Molecular coordinates
     */
    CoordMatrix * coords;
    /**
     * Atom to be minimized
     */
    const Atom * atom;

    /**
     * Step sizes for x, y and z
     */
    double steps[3];

    /**
     * Storage for coordinates
     */
    double values[3];

    /**
     * Copy simplex values to molecular atom coordinates
     */
    void copyToCoordinates(const double vals[]);

    // for use by simplex function

    int getProblemSize();
    double * getStart();
    double * getSteps();
    double getConvergenceCriteria();
    int checkConvergenceIterations();
    int maxIterations();
    double evaluate(double values[]);
    double step = 0.25;
};

/**
 * A class to minimize an entire molecule using simplex.  Note that this is normally ineffective compared to the
 * atom by atom simplex
 */
class TaffSimplexMolecule {
    friend class Simplex<TaffSimplexMolecule> ;

public:
    TaffSimplexMolecule(Taff & taff_, const int nCycles_ = 5000) :
            taff(taff_), molecule(taff_.getMolecule()), nCycles(nCycles_) {
    }

    virtual ~TaffSimplexMolecule() {
        if (steps != nullptr) {
            delete[] steps;
        }
    }

    TaffSimplexMolecule(const TaffSimplexMolecule & rhs) = delete;
    TaffSimplexMolecule & operator =(const TaffSimplexMolecule & rhs) = delete;
    TaffSimplexMolecule(TaffSimplexMolecule && rhs) = delete;
    TaffSimplexMolecule & operator =(TaffSimplexMolecule && rhs) = delete;
    /**
     * Minimize the molecule using the supplied molecular coordinates
     *
     * @param atom
     * @param coords
     * @return
     */
    double minimize(CoordMatrix & coords);

    void setStep(double step) {
          this->step = step;
      }

private:

    /**
     * Forcefield parameters for molecule
     */
    Taff & taff;

    /**
     *
     */
    const Molecule & molecule;
    /**
     * Molecular coordinates
     */
    CoordMatrix * coords;
    /**
     * Atom to be minimized
     */
    const Atom * atom;

    /**
     * Step sizes for x, y and z
     */
    double * steps = nullptr;

    double step = 0.25;

    /**
     * Storage for initial and final coordinates
     */
    Matrix<double, 3, Dynamic> atomCoords;

    // maximum number of simples cycles
    const int nCycles;

    // for use by simplex function

    int getProblemSize();
    double * getStart();
    double * getSteps();
    double getConvergenceCriteria();
    int checkConvergenceIterations();
    int maxIterations();
    double evaluate(double values[]);

};

} /* namespace GarethFF */

#endif /* SRC_FF_TAFFMINIMIZER_H_ */
