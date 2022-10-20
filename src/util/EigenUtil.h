//
// Created by gareth on 11/11/20.
//

#ifndef GAPE_C_EIGENUTIL_H
#define GAPE_C_EIGENUTIL_H

#include <cstdio>

#include <Eigen/Dense>

namespace GarethUtil {

    using namespace std;

/**
 * Prints an Eigen matrix to stdout (useful for debugging)
 *
 * @param m
 */
    void printMatrix(const Eigen::MatrixXd &m);

/**
 * Prints an Eigen vector to stdout
 *
 * @param v
 */
    void printVector(const Eigen::VectorXd &v);

}
#endif //GAPE_C_EIGENUTIL_H

