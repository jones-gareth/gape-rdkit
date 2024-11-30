//
// Created by gareth on 10/18/22.
//

#include "SuperpositionChromosome.h"

#include "RotatableBond.h"
#include "SuperpositionGa.h"
#include "util/CoordOps.h"
#include "util/LeastSquaresFit.h"
#include "util/Reporter.h"
#include <set>

namespace Gape {
    struct BestThree {
        double firstSqrDist = 1e100, secondSqrDist = 1e100, thirdSqrDist = 1e100;
        int firstNo, secondNo, thirdNo;

        void checkBest(const int index, const double sqrDist) {
            REPORT(Reporter::DEBUG) << "checking " << index << " sqrDist " << sqrDist;
            if (sqrDist < firstSqrDist) {
                thirdNo = secondNo;
                thirdSqrDist = secondSqrDist;
                secondNo = firstNo;
                secondSqrDist = firstSqrDist;
                firstNo = index;
                firstSqrDist = sqrDist;
            } else if (sqrDist < secondSqrDist) {
                thirdNo = secondNo;
                thirdSqrDist = secondSqrDist;
                secondNo = index;
                secondSqrDist = sqrDist;
            } else if (sqrDist < thirdSqrDist) {
                thirdNo = index;
                thirdSqrDist = sqrDist;
            }
        }
    };

    const double SuperpositionChromosome::passOneDistance = 2.0;
    const double SuperpositionChromosome::maxPassDistance = 5.0;

    SuperpositionChromosome::SuperpositionChromosome(const SuperpositionGa &superpositionGa) : superpositionGa(
            superpositionGa),
        binaryStringChromosome(
            superpositionGa.getSuperposition().getBinaryStringLength(),
            superpositionGa.getRng(),
            superpositionGa.getBinaryStringChromosomePolicy()),
        integerStringChromosome(
            superpositionGa.getSuperposition().getIntegerStringLength(),
            superpositionGa.getRng(),
            superpositionGa.getIntgerStringChromosomePolicy()) {
    }

    void SuperpositionChromosome::copyGene(const SuperpositionChromosome &from) {
        integerStringChromosome.copyGene(from.integerStringChromosome);
        binaryStringChromosome.copyGene(from.binaryStringChromosome);
    }

    void SuperpositionChromosome::initialize() {
        integerStringChromosome.initialize();
        binaryStringChromosome.initialize();
    }

    void SuperpositionChromosome::fitMolecules(bool remap) {
        fitted = false;

        const auto &superposition = superpositionGa.getSuperposition();
        const auto &molecules = superposition.getMolecules();
        conformerCoordinates.clear();
        conformerCoordinates.reserve(molecules.size());\

        if (operationName != OperationName::IntegerStringCrossover && operationName !=
            OperationName::IntegerStringMutate) {
            for (int moleculeNumber = 0; moleculeNumber < molecules.size(); moleculeNumber++) {
                const auto &molecule = molecules[moleculeNumber];
                auto moleculeCoordinates = std::make_shared<SuperpositionCoordinates>(
                    molecule->getReferenceCoordinates());
                if (molecule->isFixed() || molecule->isRigid()) {
                    continue;
                }
                auto pos = superposition.getBinaryEntryPoints()[moleculeNumber];
                conformerCoordinates.push_back(moleculeCoordinates);
                for (const auto &rotatableBond: molecule->getRotatableBonds()) {
                    double angle = .0;
                    switch (rotatableBond->getRotatableBondType()) {
                        case RotatableBondType::Full: {
                            const int decodedValue = binaryStringChromosome.decodeToInt(pos, 8);
                            pos += 8;
                            angle = (decodedValue * 256.0) / (2.0 * M_PI);
                        }
                        break;
                        case RotatableBondType::Flip:
                            angle = binaryStringChromosome.getValue(pos) ? M_PI : 0;
                            ++pos;
                            break;
                        default: ;
                    }
                    if (angle != .0) {
                        rotatableBond->rotateBond(angle, *moleculeCoordinates);
                    }
                }
            }
        }

        fittedCoordinates.clear();
        fittedCoordinates.reserve(molecules.size());
        int moleculeNumber = 0;
        int fittingMoleculeNumber = 0;
        const auto &fittingMoleculeCoordinates = conformerCoordinates[fittingMoleculeNumber];
        const auto fittingMolecule = superposition.getFittingMolecule();
        for (const auto &molecule: molecules) {
            const auto otherCoordinates = std::make_shared<SuperpositionCoordinates>(
                *conformerCoordinates[moleculeNumber]);
            fittedCoordinates.push_back(otherCoordinates);
            if (molecule.get() == fittingMolecule)
                continue;
            if (molecule->isFixed())
                continue;
            auto pos = superposition.getIntegerEntryPoints()[fittingMoleculeNumber++];
            ok = fitMolecule(pos, *fittingMolecule, *molecule, *fittingMoleculeCoordinates, *otherCoordinates, remap);
            if (!ok) {
                return;
            }
        }

        fitted = true;
    }

    bool SuperpositionChromosome::fitMolecule(const int start, const SuperpositionMolecule &fittingMolecule,
                                              const SuperpositionMolecule &otherMolecule,
                                              const SuperpositionCoordinates &fittingCoordinates,
                                              SuperpositionCoordinates &otherCoordinates, bool remap) {
        auto &fittingFeatures = fittingMolecule.getAllFeatures();
        int pos = start;
        std::map<std::shared_ptr<Feature>, std::shared_ptr<Feature> > featureMap;
        for (const auto &mappingFeature: fittingFeatures) {
            const int otherMapping = integerStringChromosome.getValue(pos++);
            if (otherMapping == -1)
                continue;
            const auto &otherMoleculeFeatures = otherMolecule.getFeatures();
            if (const auto &otherFeatures = otherMoleculeFeatures.find(mappingFeature->getFeatureType());
                otherFeatures != otherMoleculeFeatures.end()) {
                const auto &otherFeature = otherFeatures->second[otherMapping];
                featureMap[mappingFeature] = otherFeature;
            }
        }

        auto numberPoints = featureMap.size();
        if (numberPoints < 3) {
            REPORT(Reporter::DETAIL) << "First pass less than 3 fitting points";
            return false;
        }

        CoordMatrix fittingCoords(4, numberPoints), otherCoords(4, numberPoints);
        int pointNumber = 0;
        for (const auto &[fittingFeature, otherFeature]: featureMap) {
            auto fittingPoint = fittingFeature->getFittingPoint(fittingCoordinates);
            auto otherPoint = otherFeature->getFittingPoint(otherCoordinates);
            for (int i = 0; i < 3; i++) {
                fittingCoords(i, pointNumber) = fittingPoint[i];
                otherCoords(i, pointNumber) = otherPoint[i];
            }
            fittingCoords(3, pointNumber) = 1.0;
            otherCoords(3, pointNumber) = 1.0;
            ++pointNumber;
        }

        auto matrix = leastSquaresFit(otherCoords, fittingCoords, boost::none);
        std::set<int> closeMappings;

        if (remap) {
            auto transformedOtherCoords = matrix * otherCoords;
            REPORT(Reporter::TRACE) << "First pass transformed coordinates " << std::endl << transformedOtherCoords;
            // 2nd pass fitting
            int nMapped = 0;
            BestThree bestThree;
            for (int i = 0; i < pointNumber; i++) {
                auto c1 = transformedOtherCoords.col(i);
                auto sqrDist = sqrDistance(transformedOtherCoords.col(i), fittingCoords.col(i));
                if (sqrDist < passOneDistance * passOneDistance) {
                    closeMappings.insert(i);
                }
                bestThree.checkBest(i, sqrDist);
            }
            // If we can't find three points then use the best three points
            // within MAX_PASS_DISTANCE
            if (closeMappings.size() < 3) {
                REPORT(Reporter::DETAIL) <<
                                         "Best Three " + bestThree.firstNo << " "
                                         << bestThree.secondNo << " " << bestThree.thirdNo;
                // If the third point is too far away then the mapping has
                // failed
                if (bestThree.thirdSqrDist > maxPassDistance * maxPassDistance) {
                    REPORT(Reporter::DETAIL) << "fitMolecule 2nd pass failed";
                    return false;
                }
                closeMappings.insert(bestThree.firstNo);
                closeMappings.insert(bestThree.secondNo);
                closeMappings.insert(bestThree.thirdNo);
            }

            auto numberPoints = closeMappings.size();
            CoordMatrix fittingCoords2(4, numberPoints), otherCoords2(4, numberPoints);
            int pointNumber = 0;
            for (auto it = closeMappings.begin(); it != closeMappings.end(); ++it) {
                auto featureNumber = *it;
                auto fittingFeature = fittingFeatures[featureNumber];
                auto otherFeature = featureMap[fittingFeature];
                auto fittingPoint = fittingFeature->getFittingPoint(fittingCoordinates);
                auto otherPoint = otherFeature->getFittingPoint(otherCoordinates);
                for (int i = 0; i < 3; i++) {
                    fittingCoords(i, pointNumber) = fittingPoint[i];
                    otherCoords(i, pointNumber) = otherPoint[i];
                }
                fittingCoords(3, pointNumber) = 1.0;
                otherCoords(3, pointNumber) = 1.0;
                ++pointNumber;
            }

            matrix = leastSquaresFit(otherCoords, fittingCoords, boost::none);
        }

        otherCoordinates.transformCoordinates(matrix);
        // TODO Remap chromosome
        if (remap) {
            for (int i = 0; i < fittingFeatures.size(); i++) {
                if (closeMappings.find(i) == closeMappings.end()) {
                    auto pos = start + 1;
                    integerStringChromosome.setValue(pos, -1);
                }
                REPORT(Reporter::DEBUG) << "V " << integerStringChromosome.getValue(pos);
            }
        }

        if (Gape::Reporter::isReportingAt(Reporter::DEBUG)) {
            for (const auto &mappingFeature: fittingFeatures) {
                REPORT(Reporter::DEBUG) << "V " << integerStringChromosome.getValue(pos);
                auto pos = start + 1;
                if (auto otherFeature = featureMap.find(mappingFeature); otherFeature != featureMap.end()) {
                    REPORT(Reporter::DEBUG)
                        << mappingFeature->mappingInfo(*otherFeature->second, fittingCoordinates, otherCoordinates);
                }
            }
        }

        return true;
    }

    bool SuperpositionChromosome::isOk() {
        if (!fitted) {
            fitMolecules(true);
        }
        return ok;
    }

    double SuperpositionChromosome::score() {
        hasFitness = false;
        return getFitness();
    }

    bool SuperpositionChromosome::equals(const Gape::SuperpositionChromosome &other) const {
        if (!binaryStringChromosome.equals(other.binaryStringChromosome)) {
            return false;
        }
        return integerStringChromosome.equals(other.integerStringChromosome);
    }

    double SuperpositionChromosome::getFitness() {
        if (hasFitness)
            return fitness;

        assert(ok && fitted);

        fitness = .0;
        const auto &settings = superpositionGa.getSuperposition().settings;
        // Conformational Energies
        if (settings.getGapeParameters().conformationalWeight > .0) {
            calculateConformationalEnergy();
        }

        // Volume Overlay
        if (settings.getGapeParameters().volumeWeight > 0) {
            calculateVolumeIntegral();
        }


        // All the features now:

        // Initialize pharmacophore information
        featureMapping.clear();
        const auto &molecules = superpositionGa.getSuperposition().getMolecules();
        for (int moleculeNumber = 0; moleculeNumber < molecules.size(); moleculeNumber++) {
            const auto & molecule = molecules[moleculeNumber];
            const auto superpositionCoordinates = conformerCoordinates[moleculeNumber];
            for (const auto feature : molecule->getAllFeatures()) {
                const auto ptr = feature.get();
                const auto point = feature->getFittingPoint(*superpositionCoordinates);
                auto featureInformation = std::make_shared<FeatureInformation>(ptr, point);
                featureMapping[feature.get()] = featureInformation;
            }
        }

        // For donors, acceptors and aromatic rings given the standard
        // configuration of compareAll not set and singleMatchOnly set we'll end
        // up calling singleMatchFeatureOverlay for each of these.

        // Donor Hydrogens

        if (settings.getGapeParameters().donorHydrogenWeight > .0) {

        }
        /*
        if (superposition.getDonorHydrogenWt() > 0)
            donorHydrogenScore();

        // Acceptor Atom
        if (superposition.getAcceptorAtomWt() > 0)
            acceptorAtomScore();

        // Aromatic Ring
        if (superposition.getAromaticRingWt() > 0)
            aromaticRingScore();

        // Constraints (including breaking large cycles)
        if (superposition.getConstraintWt() > 0)
            constraintScore();

        fitness = superposition.getDonorHydrogenWt() * donorHydrogenScore
                  + superposition.getAcceptorAtomWt() * acceptorAtomScore
                  + superposition.getAromaticRingWt() * aromaticRingScore
                  + superposition.getVolumeWt() * volumeIntegral
                  - superposition.getConfWt() * conformationalEnergy
                  - superposition.getConstraintWt() * constraintScore;

        // UserFeatures
        for (int i = 0; i < superposition.getnUserFeatureTypes(); i++) {
            FeatureType featureType = Feature.userFeatureType(i);
            userFeatureScore(featureType);
            UserFeatureSet featureSet = superposition.getUserFeatureSet(featureType);
            fitness += featureSet.getFeatureWeight() * userFeatureScores[i];
        }

        logger.debug("Fitness " + fitness);

         */
        hasFitness = true;
        return fitness;
    }

    double SuperpositionChromosome::rebuild() {
        hasFitness = false;
        fitMolecules(false);
        return getFitness();
    }

    double SuperpositionChromosome::calculateConformationalEnergy() {
        double conformationalEnergy = .0;
        const auto &molecules = superpositionGa.getSuperposition().getMolecules();
        conformationalEnergies.clear();
        conformationalEnergies.reserve(molecules.size());
        for (int i = 0; i < molecules.size(); i++) {
            const auto &molecule = molecules[i];
            const auto &conformer = conformerCoordinates[i]->getConformer();
            const auto moleculeConformerEnergy = molecule->calculateConformationalEnergy(conformer);
            conformationalEnergies.push_back(moleculeConformerEnergy);
            conformationalEnergy += moleculeConformerEnergy;
            REPORT(Reporter::DEBUG) << "conformational energy for molecule " << i << " " << moleculeConformerEnergy;
        }
        conformationalEnergy = conformationalEnergy / molecules.size();
        return conformationalEnergy;
    }

    double SuperpositionChromosome::calculateVolumeIntegral() {
    }
}
