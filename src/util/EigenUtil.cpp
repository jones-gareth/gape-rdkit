//
// Created by gareth on 11/11/20.
//

#include "EigenUtil.h"
#include <iostream>

namespace Gape {

    using namespace std;

    void printMatrix(const Eigen::MatrixXd &m) {
        cout << "Matrix " << endl << m;
    }

    void printVector(const Eigen::VectorXd &v) {
        cout << "Vector " << endl << v;

    }

}