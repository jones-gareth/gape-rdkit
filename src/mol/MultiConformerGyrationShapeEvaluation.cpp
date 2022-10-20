/*
 * MultiConformerGyrationShapeEvaluation.cpp
 *
 *  Created on: Apr 12, 2016
 *      Author: gjones
 */

#include "MultiConformerGyrationShapeEvaluation.h"
#include "../util/ShapeMetrics.h"
#include "../util/Reporter.h"
#include "../util/DetermineNearestNeighbourLists.h"

namespace GarethMol {

void MultiConformerGyrationShapeEvalution::buildDescriptorValues() {
    auto nDescriptors = 3;
    auto nConformers = molecule.nConformers();
    descriptorValues = make_unique<Array2D<double>>(nConformers, nDescriptors);
    auto nHeavyAtoms = molecule.nHeavyAtoms();

    for (auto i = 0ul; i < nConformers; i++) {
        // calculate shape descriptors on heavy atoms only
        auto heavyCoords = CoordMatrix(4, nHeavyAtoms);
        auto allCoords = molecule.getConformer(i).getCoordinates();
        auto heavyAtomNo = 0;
        for (auto i = 0ul; i < molecule.nAtoms(); i++) {
            if (AtomType::isHeavy(molecule.getAtom(i).getAtomTypeId())) {
                heavyCoords.col(heavyAtomNo) = allCoords.col(i);
                heavyAtomNo++;
            }
        }
        assert(heavyAtomNo == nHeavyAtoms);
        ShapeMetrics metrics(heavyCoords);
        if (false) {
            auto values = metrics.getTensorSingularValues();
            descriptorValues->set(i, 0, values[0]);
            descriptorValues->set(i, 1, values[1]);
            descriptorValues->set(i, 2, values[2]);
        } else {
            descriptorValues->set(i, 0, metrics.getRadiusOfGyration());
            descriptorValues->set(i, 1, metrics.getAsphericity());
            descriptorValues->set(i, 2, metrics.getAcylindricity());
            //     descriptorValues->set(i, 3, metrics.getAnisotrophy());
        }
    }

    descriptorValues->normalizeColumns();

    for (auto i = 0ul; i < nConformers; i++) {
        for (auto j = 0ul; j < nConformers; j++) {
            auto sqrDistance = descriptorValues->rowSqrDistance(i, j);
            auto distance = sqrt(sqrDistance);
            distanceMatrix.set(i, j, distance);
            distanceMatrix.set(j, i, distance);
            REPORT(Reporter::DEBUG) << "Normalized shape distance between "
                    << to_string(i + 1) << " and " << to_string(j + 1) << " is "
                    << distance;
        }
    }
}

const std::vector<std::vector<size_t> > MultiConformerGyrationShapeEvalution::determineNearestNeighbourLists() {
    nearestNeighbours =
            DetermineNearestNeighbourLists<MultiConformerGyrationShapeEvalution>::determineNearestNeighbourLists(
                    *this);
    return nearestNeighbours;
}

} /* namespace GarethMol */
