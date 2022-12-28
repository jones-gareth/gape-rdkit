//
// Created by jones on 12/23/2022.
//

#include "TransformOps.h"


#include <Geometry/point.h>

 namespace Gape {

     using namespace RDKit;

    void determineRotation(const RDGeom::Point3D &point1, const RDGeom::Point3D &point2, double angle, RDGeom::Transform3D &transform) {
        RDGeom::Point3D diff = point1 - point2;
        diff.normalize();
        RDGeom::Transform3D trans, xRot, yRot, rot, invTrans, invXRot, invYRot;
        trans.SetTranslation(-point1);
        invTrans.SetTranslation(point1);
        rot.SetRotation(angle, diff);
        transform.assign(invTrans*rot*trans);
    }

}