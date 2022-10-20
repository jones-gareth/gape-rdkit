/*
 * CheckTypesTest.cxx
 *
 *	Test file for automatic assignment of atom types
 *
 *  Needs some more tests
 *
 */

#include <gtest/gtest.h>

#include "../util/CoordOps.h"
#include "../util/Reporter.h"

namespace {

using namespace std;
using namespace GarethUtil;
using namespace GarethMol;
using namespace ::testing;

/**
 * Test class for some geometric transforms
 */
class GeometryTest: public Test {
protected:

	GeometryTest() {
		Reporter::setMinReportingLevel(Reporter::NORMAL);
	}
};

/**
 * Test the transformation for moving point1 to the origin and point2 to the -ve z axis is working fine
 * @param point1
 * @param point2
 */
void testTransform(const CoordVector &point1, const CoordVector & point2) {
	Transform<double, 3, Affine> in, out;
	transformToNegativeZaxis(point1, point2, in, out);

	auto newPoint = in * point2;
	auto checkPoint = out * newPoint;

	REPORT(Reporter::DEBUG) << "translated " << point2 << " to " << newPoint
			<< " and back to " << checkPoint;

	// static casts to prevent eclipse cdt ide error

	ASSERT_TRUE(equals(static_cast<double>(newPoint(0)), 0.0, 1e-10));
	ASSERT_TRUE(equals(static_cast<double>(newPoint(1)), 0.0, 1e-10));
	ASSERT_TRUE(static_cast<double>(newPoint(2)) < .0);
	for (int i = 0; i < 3; i++) {
		ASSERT_TRUE(
				equals(static_cast<double>(checkPoint(i)),
						static_cast<double>(point2(i)), 1e-10));
	}

	newPoint = in * point1;
	checkPoint = out * newPoint;
	ASSERT_TRUE(equals(static_cast<double>(newPoint(0)), 0.0, 1e-10));
	ASSERT_TRUE(equals(static_cast<double>(newPoint(1)), 0.0, 1e-10));
	ASSERT_TRUE(equals(static_cast<double>(newPoint(2)), 0.0, 1e-10));
	for (int i = 0; i < 3; i++) {
		ASSERT_TRUE(
				equals(static_cast<double>(checkPoint(i)),
						static_cast<double>(point1(i)), 1e-10));
	}
}

/**
 * Check set up trans actually moves coordinates to the -ve z axis
 */
TEST_F(GeometryTest, testSetUpTrans) {
	testTransform(CoordVector(0.0, 0.0, 0.0, 1.0),
			CoordVector(0.0, 0.0, 1.0, 1.0));
	testTransform(CoordVector(0.0, 0.0, 0.0, 1.0),
			CoordVector(0.0, 1.0, 0.0, 1.0));
	testTransform(CoordVector(0.0, 0.0, 0.0, 1.0),
			CoordVector(1.0, 0.0, 0.0, 1.0));
	testTransform(CoordVector(0.0, 0.0, 0.0, 1.0),
			CoordVector(1.0, 1.0, 1.0, 1.0));
	testTransform(CoordVector(0.0, 0.0, 0.0, 1.0),
			CoordVector(-1.0, 0.0, 0.0, 1.0));
	testTransform(CoordVector(0.0, 0.0, 0.0, 1.0),
			CoordVector(0.0, -1.0, 0.0, 1.0));
	testTransform(CoordVector(0.0, 0.0, 0.0, 1.0),
			CoordVector(0.0, 0.0, -1.0, 1.0));

	testTransform(CoordVector(1.0, 1.0, 1.0, 1.0),
			CoordVector(1.0, 1.0, 2.0, 1.0));
	testTransform(CoordVector(1.0, 1.0, 1.0, 1.0),
			CoordVector(-1.0, -1.0, -1.0, 1.0));

	for (double i = -1; i < 2; i = i + 0.25) {
		for (double j = -1; j < 2; j = j + 0.25) {
			for (double k = -1; k < 2; k = k + 0.25) {
				if (i == 0 && j == 0 && k == 0)
					continue;
				testTransform(CoordVector(0.0, 0.0, 0.0, 1.0),
						CoordVector(i, j, k, 1.0));
			}
		}
	}

}

}
