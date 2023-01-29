/*
 * LeastSquaresFit.cpp
 *
 *  Created on: Sep 9, 2015
 *      Author: Gareth Jones
 */

#include <cassert>

#include <Eigen/SVD>

#include "LeastSquaresFit.h"
#include "Reporter.h"

namespace Gape {

using namespace std;
using namespace Eigen;

// set true to use Eigen Affine transforms
static const bool USE_AFFINE = false;

const Transform<double, 3, Affine> leastSquaresFit(const MatrixXd & xCoords,
		const MatrixXd & yCoords,
		const boost::optional<const VectorXd &> weights) {

	auto nPoints = xCoords.cols();
	auto hasWeights = weights.is_initialized();

	assert(yCoords.cols() == nPoints);
	assert(xCoords.rows() == 3 || xCoords.rows() == 4);
	assert(yCoords.rows() == 3 || yCoords.rows() == 4);
	if (hasWeights) {
		assert(weights->size() == nPoints);
	}

	// create weights matrix (or identity if no weights matrix)
	MatrixXd weightsMatrix(nPoints, nPoints);
	if (hasWeights) {
		weightsMatrix = weights->asDiagonal();
	} else {
		weightsMatrix.setIdentity();
	}
	REPORT(Reporter::TRACE) << "Weights matrix " << endl << weightsMatrix;

	// create nPoints by 3 coordinate matrices
	MatrixXd x = xCoords.block(0, 0, 3, nPoints);
	MatrixXd y = yCoords.block(0, 0, 3, nPoints);
	REPORT(Reporter::TRACE) << "X " << endl << x;
	REPORT(Reporter::TRACE) << "Y " << endl << y;

	// create weighted matrices
	MatrixXd weightedXMatrix = x * weightsMatrix;
	MatrixXd weightedYMatrix = y * weightsMatrix;
	REPORT(Reporter::TRACE) << "weighted X " << endl << weightedXMatrix;
	REPORT(Reporter::TRACE) << "weighted Y " << endl << weightedYMatrix;
	REPORT(Reporter::TRACE) << "X " << endl << x << endl << "Y " << endl << y;

	// determine X and Y weighted centers
	auto weightSum = weightsMatrix.sum();
	REPORT(Reporter::TRACE) << "weighted sum " << weightSum;
	auto xCenter = (weightedXMatrix.rowwise().sum()) / weightSum;
	auto yCenter = (weightedYMatrix.rowwise().sum()) / weightSum;
	REPORT(Reporter::TRACE) << "X center" << endl << xCenter << endl
			<< "Y center" << endl << yCenter;

	// center input data
	for (auto i = 0; i < nPoints; i++) {
		x.col(i) -= xCenter;
		y.col(i) -= yCenter;
	}
	REPORT(Reporter::TRACE) << "X centered" << endl << x << endl;
	REPORT(Reporter::TRACE) << "Y centered" << endl << x << endl;

	MatrixXd XYPrime = x * weightsMatrix * y.transpose();
	REPORT(Reporter::TRACE) << "Product matrix prior to SVD " << endl
			<< XYPrime;

	JacobiSVD<Matrix3d> svd(XYPrime, ComputeFullU | ComputeFullV);
	Matrix3d uPrime = svd.matrixU().transpose();
	REPORT(Reporter::TRACE) << "U post SVD " << endl << svd.matrixU();
	REPORT(Reporter::TRACE) << "V post SVD " << endl << svd.matrixV();
	REPORT(Reporter::TRACE) << "Singular values " << endl
			<< svd.singularValues();

	Matrix3d rotMatrix = svd.matrixV() * uPrime;
	auto det = rotMatrix.determinant();
	REPORT(Reporter::TRACE) << "Determinant is " << det;

	if (det < .0) {
		// If the determinant is 1 then the matrix is only
		// rotation. Otherwise we look through the singular values
		// (they are not ordered) and change the sign of the row
		// in U' that corresponds to the smallest singular value.

		// in Eigen SVD the singular values are always in descreasing order
		uPrime.row(2) *= -1;
		rotMatrix = svd.matrixV() * uPrime;
		det = rotMatrix.determinant();
		REPORT(Reporter::TRACE) << "Fixed determinant is " << det;
		assert(det > 0);
	}

	Matrix4d rotate;
	rotate.setZero();
	rotate(3, 3) = 1.0;
	rotate.block(0, 0, 3, 3) = rotMatrix;
	REPORT(Reporter::TRACE) << "Rotation matrix " << endl << rotate;

	if (USE_AFFINE) {
		// using Eigen affine transformations reduces a lot of the biolerplate and is clearer
		Transform<double, 3, Affine> rotateTrans, inTrans, outTrans;
		rotateTrans.matrix() = rotate;
		inTrans = Translation<double, 3>(-xCenter);
		REPORT(Reporter::TRACE) << "Affine in translation matrix " << endl
				<< inTrans.matrix();
		outTrans = Translation<double, 3>(yCenter);
		REPORT(Reporter::TRACE) << "Affine out translation matrix " << endl
				<< outTrans.matrix();
		Transform<double, 3, Affine> trans = outTrans * rotateTrans * inTrans;
		return trans;
	} else {
		// otherwise this is the old way of doing things
		Matrix4d in;
		in.setIdentity();
		in.block(0, 3, 3, 1) = -xCenter;
		REPORT(Reporter::TRACE) << "In translation matrix " << endl << in;

		Matrix4d out;
		out.setIdentity();
		out.block(0, 3, 3, 1) = yCenter;
		REPORT(Reporter::TRACE) << "Out translation matrix " << endl << out;

		Matrix4d matrix = out * rotate * in;
		REPORT(Reporter::TRACE) << "Final product matrix " << endl << matrix;

		// if we return to this probably just want to return the matrix
		Transform<double, 3, Affine> transform;
		transform.matrix() = matrix;
		return transform;
	}
}

/**
 * Print the matrix as a string.
 *
 * @param matrix
 * @return
 */
string matrixToString(MatrixXd matrix) {
	return printToString(matrix);
}

} /* namespace Gape */
