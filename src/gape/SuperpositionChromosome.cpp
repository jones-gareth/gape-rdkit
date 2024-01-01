//
// Created by gareth on 10/18/22.
//

#include "SuperpositionChromosome.h"

#include "RotatableBond.h"
#include "SuperpositionGa.h"
#include "util/CoordOps.h"
#include "util/LeastSquaresFit.h"
#include "util/Reporter.h"

namespace Gape {
    SuperpositionChromosome::SuperpositionChromosome(const SuperpositionGa& superpositionGa): superpositionGa(
            superpositionGa),
        binaryStringChromosome(superpositionGa.getSuperposition().getBinaryStringLength(), superpositionGa.getRng(),
                               superpositionGa.getBinaryStringChromosomePolicy()),
        integerStringChromosome(superpositionGa.getSuperposition().getIntegerStringLength(), superpositionGa.getRng(),
                                superpositionGa.getIntgerStringChromosomePolicy()) {
    }

    void SuperpositionChromosome::copyGene(const SuperpositionChromosome& from) {
        integerStringChromosome.copyGene(from.integerStringChromosome);
        binaryStringChromosome.copyGene(from.binaryStringChromosome);
    }

    void SuperpositionChromosome::initialize() {
        integerStringChromosome.initialize();
        binaryStringChromosome.initialize();
    }

    double SuperpositionChromosome::score() {
    }

    void SuperpositionChromosome::fitMolecules(bool remap) {
        fitted = false;

        const auto& superposition = superpositionGa.getSuperposition();
        const auto& molecules = superposition.getMolecules();
        conformerCoordinates.clear();
        conformerCoordinates.reserve(molecules.size());\

        if (operationName != OperationName::IntegerStringCrossover && operationName !=
            OperationName::IntegerStringMutate) {
            for (int moleculeNumber = 0; moleculeNumber < molecules.size(); moleculeNumber++) {
                const auto& molecule = molecules[moleculeNumber];
                auto pos = superposition.getBinaryEntryPoints()[moleculeNumber];
                auto moleculeCoordinates = std::make_shared<SuperpositionCoordinates>(
                    molecule->getReferenceCoordinates());
                conformerCoordinates.push_back(moleculeCoordinates);
                for (const auto& rotatableBond: molecule->getRotatableBonds()) {
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
                        default:
                            ;
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
        const auto& fittingMoleculeCoordinates = conformerCoordinates[fittingMoleculeNumber];
        for (const auto& molecule: molecules) {
            fittedCoordinates.push_back(conformerCoordinates[moleculeNumber++]);
            if (molecule.get() == superposition.getFittingMolecule())
                continue;
            if (molecule->isFixed())
                continue;
            auto pos = superposition.getIntegerEntryPoints()[fittingMoleculeNumber++];
        }
    }

    bool SuperpositionChromosome::fitMolecule(const int start, const SuperpositionMolecule& fittingMolecule,
                                              const SuperpositionMolecule& otherMolecule,
                                              const SuperpositionCoordinates& fittingCoordinates,
                                              SuperpositionCoordinates otherCoordinates, bool remap) {
        auto& fittingFeatures = fittingMolecule.getAllFeatures();
        int pos = start;
        std::map<std::shared_ptr<Feature>, std::shared_ptr<Feature>> featureMap;
        for (const auto& mappingFeature: fittingFeatures) {
            const int otherMapping = integerStringChromosome.getValue(pos++);
            if (otherMapping == -1)
                continue;
            const auto& otherMoleculeFeatures = otherMolecule.getFeatures();
            if (const auto& otherFeatures = otherMoleculeFeatures.find(mappingFeature->getFeatureType());
                otherFeatures != otherMoleculeFeatures.end()) {
                const auto& otherFeature = otherFeatures->second[otherMapping];
                featureMap[mappingFeature] = otherFeature;
            }
        }

        auto numberPoints = featureMap.size();
        if (numberPoints < 3) {
            REPORT(Reporter::DEBUG) << "First pass less than 3 fitting points";
            return false;
        }

        CoordMatrix fittingCoords(4, numberPoints), otherCoords(4, numberPoints);
        int pointNumber = 0;
        for (const auto& [fittingFeature, otherFeature]: featureMap) {
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
        auto transformedOtherCoords = matrix * otherCoords;
        REPORT(Reporter::TRACE) << "First pass transformed corordinates " << std::endl << transformedOtherCoords;


        return true;
    }
}
