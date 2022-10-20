/*
 * ShapeMetricsMoleculeTest.cpp
 *
 *  Simple tests for calculating shape metrics
 *
 *  Created on: Sep 10, 2015
 *      Author: gjones
 */

#include <gtest/gtest.h>
#include "../mol/MulticonformerMolecule.h"
#include "../mol/MultiConformerGyrationShapeEvaluation.h"
#include "../util/CoordOps.h"
#include "../util/ShapeMetrics.h"
#include "../util/Reporter.h"
#include "../util/Util.h"

namespace {

using namespace std;
using namespace GarethUtil;
using namespace ::testing;

class ShapeMetricsMoleculeTest: public Test {

protected:

    static void SetUpTestCase() {
        Reporter::setMinReportingLevel(Reporter::DEBUG);

        //Reporter::setFileStream(out);
    }

    static void TearDownTestCase() {
        out.close();
    }

    static ofstream out;
};

ofstream ShapeMetricsMoleculeTest::out { "ShapeMetricsMoleculeTest.log" };

/**
 * Simple test- just confirms that we get the same results for the radius of
 * gyration using different methods
 *
 * @param inFile
 */
void testLigandConformerFile(const string & inFile) {

    REPORT(Reporter::INFO) << "loading file " << inFile;

    auto molecule = MulticonformerMolecule::readAndInitializeMolecule(inFile);

    for (const auto & conformer : molecule->getConformers()) {

        double radiusOfGyration = ::radiusOfGyration(
                conformer->getCoordinates());
        ShapeMetrics shapeMetrics(conformer->getCoordinates());

        REPORT(Reporter::NORMAL) << "Radius of gyration " << radiusOfGyration;
        REPORT(Reporter::NORMAL) << "Radius of gyration from gyration tensor "
                << shapeMetrics.getRadiusOfGyration();
        REPORT(Reporter::NORMAL) << "Asphericity   "
                << shapeMetrics.getAsphericity();
        REPORT(Reporter::NORMAL) << "Acylindricity "
                << shapeMetrics.getAcylindricity();
        REPORT(Reporter::NORMAL) << "Anisotrophy   "
                << shapeMetrics.getAnisotrophy();

        ASSERT_TRUE(
                equals(radiusOfGyration, shapeMetrics.getRadiusOfGyration(),
                        1e-6));
    }


    MultiConformerGyrationShapeEvalution eval(*molecule);
    eval.determineNearestNeighbourLists();

}

/**
 * Test of ligand loading
 */
TEST_F(ShapeMetricsMoleculeTest, Test1g9v) {
    testLigandConformerFile("../resources/1g9v_omega_20.sdf");
}

}
