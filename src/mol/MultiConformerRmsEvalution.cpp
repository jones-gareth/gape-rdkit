/*
 * MultiConformerRmsEvalution.cpp
 *
 *  Created on: Sep 15, 2015
 *      Author: gjones
 */

#include <vector>
#include <limits>

#include "MultiConformerRmsEvalution.h"
#include "MolComparer.h"
#include "MolecularRms.h"
#include "../util/DetermineNearestNeighbourLists.h"

namespace GarethMol {

    using namespace GarethUtil;
    using namespace std;

    const Array2D<double> &MultiConformerRmsEvalution::evaluateRms(bool matchTypes) {

        auto molComparer = MolComparer::createIsomorphismComparer(molecule);
        molComparer->setHeavyAtomOnly(true);
        molComparer->setSubgraph(false);
        molComparer->setMatchTypes(matchTypes);

        auto callbackFunction =
                [this](const std::vector<size_t> &queryIdsToTargetIds) {
                    addIsomorphism(queryIdsToTargetIds);
                };
        molComparer->setCallbackFunction(callbackFunction);
        molComparer->compare();

        const auto nConformers = molecule.nConformers();
        for (auto i = 0ul; i < nConformers; i++) {

            REPORT(Reporter::DEBUG)
                << "Determining conformer rms for conformer number "
                << to_string(i + 1);
            rmsValues(i, i) = 0.0;

            for (auto j = i + 1; j < nConformers; j++) {

                const auto &conformer1 = molecule.getConformers().at(i);
                const auto &conformer2 = molecule.getConformers().at(j);
                const auto &coords1 = conformer1->getCoordinates();
                const auto &coords2 = conformer2->getCoordinates();

                double minRms = numeric_limits<double>::max();
                for (auto mapping : isomorphisms) {

                    double rms = MolecularRms::leastSquaresAtomsFit(mapping,
                                                                    molComparer->getQueryAtoms(),
                                                                    molComparer->getTargetAtoms(), coords1, coords2);
                    if (rms < minRms)
                        minRms = rms;
                }

                rmsValues(i, j) = minRms;
                rmsValues(j, i) = minRms;
            }
        }

        return rmsValues;
    }

    void MultiConformerRmsEvalution::addIsomorphism(
            const std::vector<size_t> &queryIdsToTargetIds) {
        // add isomorphism to list
        isomorphisms.push_back(queryIdsToTargetIds);
    }

    const std::vector<std::vector<size_t>> &MultiConformerRmsEvalution::determineNearestNeighbourLists() {
        nearestNeighbours =
                DetermineNearestNeighbourLists<MultiConformerRmsEvalution>::determineNearestNeighbourLists(
                        *this);
        return nearestNeighbours;
    }

    bool MultiConformerRmsEvalution::checkNeighbourLists() const {
        return DetermineNearestNeighbourLists<MultiConformerRmsEvalution>::checkNeighbourLists(*this, nearestNeighbours);
    }

} /* namespace GarethUtil */
