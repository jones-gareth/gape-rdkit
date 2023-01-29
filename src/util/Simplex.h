/*
 * Simplex.h
 *
 *  Created on: May 6, 2016
 *      Author: Gareth Jones
 */

#ifndef SRC_UTIL_SIMPLEX_H_
#define SRC_UTIL_SIMPLEX_H_

#include <stdexcept>
#include <functional>

namespace Gape {


using namespace std;

void nelmin(function<double(double [])> fn, int n, double start[], double xmin[],
        double *ynewlo, double reqmin, double step[], int konvge, int kcount,
        int *icount, int *numres, int *ifault);


/**
 * A class to perform simplex optimization
 *
 * This is a wrapper round ASA047.
 *
 * ASA047 is Applied Statistics Algorithm 47. Source code for many Applied Statistics Algorithms
 * is available through STATLIB. Distributed under LGPL license
 *
 * http://people.sc.fsu.edu/~jburkardt/cpp_src/asa047/asa047.html
 *
 * The SimplexInfo class needs to provide the following:
 *
 * int getProblemSize()  the number of variables
 * double * getStart() starting position of variables
 * double * getSteps() stepsize for each variable
 * double getConvergenceCriteria() variance in target function for convergence
 * int checkConvergenceIterations() interval to check convergence
 * int maxIterations() maximum number of iterations
 * double info.evaluate(double values[]) get score for values
 */
template<typename SimplexInfo>
class Simplex {
public:

    enum class Error {
        NO_ERROR, ILLEGAL_VALUE, NO_CONVERGENCE
    };

    Simplex(SimplexInfo & info_) :
            info(info_) {
    }

    virtual ~Simplex() {
        if (minValues != nullptr)
            delete[] minValues;
    }

    Simplex(const Simplex & rhs) = delete;
    Simplex & operator =(const Simplex & rhs) = delete;
    Simplex(Simplex && rhs) = delete;
    Simplex & operator =(Simplex && rhs) = delete;

    /**
     * Performs the optimization.  Returns minimum score.  Optimal values are
     * placed in the staring array
     *
     * @return
     */
    double minimize() {
        int size = info.getProblemSize();
        double * start = info.getStart();
        double * steps = info.getSteps();
        double reqMin = info.getConvergenceCriteria();
        int checkConvergenceIterations = info.checkConvergenceIterations();
        int maxIterations = info.maxIterations();

        int fault;
        minValues = new double[size];

        auto func = [this] (double values[]) -> double {
            return info.evaluate(values);
        };

        nelmin(func, size, start, minValues, &minScore, reqMin, steps,
                checkConvergenceIterations, maxIterations, &nIterations,
                &numRestarts, &fault);

        switch (fault) {
        case 0:
            error = Error::NO_ERROR;
            break;
        case 1:
            error = Error::ILLEGAL_VALUE;
            break;
        case 2:
            error = Error::NO_CONVERGENCE;
            break;
        default:
            throw runtime_error("Unknown simplex fault value");
        }

        // put minimized values back into starting values
        for (auto i =0; i<size; i++) {
            start[i]=minValues[i];
        }

        return minScore;
    }

    // getters, available post optimization

    Error getError() const {
        return error;
    }

    double getMinScore() const {
        return minScore;
    }

    double* getMinValues() const {
        return minValues;
    }

    int getIterations() const {
        return nIterations;
    }

    int getNumRestarts() const {
        return numRestarts;
    }

private:
    SimplexInfo & info;
    int nIterations, numRestarts;
    double * minValues = nullptr;
    double minScore;
    Error error = Error::NO_ERROR;
};

} /* namespace Gape */

#endif /* SRC_UTIL_SIMPLEX_H_ */
