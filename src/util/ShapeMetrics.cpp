/*
 * ShapeMetrics.cpp
 *
 *  Created on: Apr 12, 2016
 *      Author: gjones
 */

#include "Reporter.h"
#include "ShapeMetrics.h"
#include <Eigen/SVD>

namespace Gape {

const Matrix3d ShapeMetrics::determineGyrationTensor(
        const CoordMatrix & coords) {
    Matrix3d tensor;
    auto size = coords.cols();
    auto dSize = static_cast<double>(size);
    double norm = 1.0 / (2.0 * dSize * dSize);

    for (auto m = 0; m < 3; m++) {
        for (auto n = 0; n < 3; n++) {

            auto sum = .0;
            for (auto i = 0; i < size; i++) {
                for (auto j = 0; j < size; j++) {
                    sum += (coords(m, i) - coords(m, j))
                            * (coords(n, i) - coords(n, j));
                }
            }
            tensor(m, n) = norm * sum;
        }
    }

    REPORT(Reporter::DEBUG) << "Gyration tensor is " << endl << tensor;
    return tensor;
}

void ShapeMetrics::calculateGyrationTensorShapeDescriptors() {
    tensor = determineGyrationTensor(coords);
    //JacobiSVD<Matrix3d> svd(tensor, ComputeThinU | ComputeThinV);
    // tensorSingularValues = svd.singularValues().reverse();

    // The gyration tensor is symmetric, so we can use the self-adjoint solver to diagonalize
    SelfAdjointEigenSolver<Matrix3d> solver(tensor);
    if (solver.info() != Success) {
        string message = "Failed to determine eigenvalues of gyration tensor";
        REPORT(Reporter::FATAL) << message;
        throw runtime_error(message);
    }
    tensorSingularValues = solver.eigenvalues();
    REPORT(Reporter::DEBUG) << "Singular values " << tensorSingularValues;
    assert(tensorSingularValues[0] >= .0);
    assert(tensorSingularValues[1] >= .0);
    assert(tensorSingularValues[2] >= .0);
    assert(tensorSingularValues[0] <= tensorSingularValues[1]);
    assert(tensorSingularValues[1] <= tensorSingularValues[2]);

    radiusOfGyration = tensorSingularValues.sum();

    asphericity = 1.5 * tensorSingularValues[2] - radiusOfGyration / 2.0;
    assert(asphericity >= 0);

    acylindricity = tensorSingularValues[1] - tensorSingularValues[0];
    assert(acylindricity >= 0);

    anisotrophy = (asphericity * asphericity
            + (3.0 * acylindricity * acylindricity) / 4.0)
            / (radiusOfGyration * radiusOfGyration);
    assert(anisotrophy >= .0);
    assert(anisotrophy <= 1.0);
}

} /* namespace Gape */
