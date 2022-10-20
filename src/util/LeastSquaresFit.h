/*
 * LeastSquaresFit.h
 *
 *  Created on: Sep 9, 2015
 *      Author: gjones
 */

#ifndef SRC_UTIL_LEASTSQUARESFIT_H_
#define SRC_UTIL_LEASTSQUARESFIT_H_

#include <boost/optional.hpp>
#include <Eigen/Dense>

namespace GarethUtil {

using namespace std;
using namespace Eigen;

const Transform<double, 3, Affine>  leastSquaresFit(const MatrixXd & xCoords, const MatrixXd & yCoords,
		const boost::optional<const VectorXd &> weights);

string matrixToString(MatrixXd matrix);

} /* namespace GarethUtil */

#endif /* SRC_UTIL_LEASTSQUARESFIT_H_ */
