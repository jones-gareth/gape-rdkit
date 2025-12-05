//
// Created by Gareth Jones on 2/23/2025.
//

#define CATCH_CONFIG_MAIN

#include <Geometry/point.h>

#include "catch2/catch.hpp"
#include <Geometry/Transform3D.h>
#include <util/Reporter.h>
#include <util/TransformOps.h>
#include <util/Util.h>

using namespace Gape;

/**
 * Test the transformation for moving origin to the origin and atom1 to the -ve z axis is working fine
 * @param point1
 * @param point2
 */
void testTransform(const RDGeom::Point3D &origin, const RDGeom::Point3D &atom1) {
    // Set up the transformation that moves origin to the origin
    // and puts atom1 on the z-axis
    RDGeom::Transform3D transformIn, transformOut, rotation;
    setUpZAxisTransformation(origin, atom1, transformIn, transformOut);

    RDGeom::Point3D movedOrigin(origin);
    transformIn.TransformPoint(movedOrigin);
    RDGeom::Point3D checkOrigin(movedOrigin);
    transformOut.TransformPoint(checkOrigin);

    RDGeom::Point3D movedAtom1(atom1);
    transformIn.TransformPoint(movedAtom1);
    RDGeom::Point3D checkAtom1(movedAtom1);
    transformOut.TransformPoint(checkAtom1);

    REPORT(Reporter::DEBUG) << "translated " << atom1 << " to " << movedAtom1
            << " and back to " << checkAtom1;

    for (int i = 0; i < 3; i++) {
        CHECK(equals(movedOrigin[i], 0.0, 1e-10));
        CHECK(equals(checkOrigin[i], origin[i], 1e-10));
    }

    CHECK(equals(movedAtom1[0], 0.0, 1e-10));
    CHECK(equals(movedAtom1[1], 0.0, 1e-10));
    CHECK(movedAtom1[2] < 0);
    for (int i = 0; i < 3; i++) {
        CHECK(equals(checkAtom1[i], atom1[i], 1e-10));
    }
}

TEST_CASE("Transform to origin and -ve Z axis", "[GapeGeometry]") {

    Reporter::setMinReportingLevel(Reporter::TRACE);
    testTransform(RDGeom::Point3D(0, 0, 0), RDGeom::Point3D(1.0, 1.0, 1.0));
    testTransform(RDGeom::Point3D(0, 0, 0), RDGeom::Point3D(0, 1.0, 0));
    testTransform(RDGeom::Point3D(0, 0, 0), RDGeom::Point3D(0, 0, 1.0));
    testTransform(RDGeom::Point3D(0, 0, 0), RDGeom::Point3D(1.0, 0, 0));
    testTransform(RDGeom::Point3D(0, 0, 0), RDGeom::Point3D(0, 0, -1.0));
    testTransform(RDGeom::Point3D(0, 0, 0), RDGeom::Point3D(0, -1.0, 0));
    testTransform(RDGeom::Point3D(0, 0, 0), RDGeom::Point3D(-1.0, 0, 0));

    testTransform(RDGeom::Point3D(1.0, 1.0, 1.0), RDGeom::Point3D(1.0, 1.0, 2.0));
    testTransform(RDGeom::Point3D(1.0, 1.0, 1.0), RDGeom::Point3D(-1.0, -1.0, -1.0));

}
