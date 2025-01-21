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
#include <GraphMol/FileParsers/MolWriters.h>

#include "FeatureOverlay.h"

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
            superpositionGa.getIntegerStringChromosomePolicy()) {
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
            BestThree bestThree;
            for (int i = 0; i < pointNumber; i++) {
                // auto c1 = transformedOtherCoords.col(i);
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

            auto _numberPoints = closeMappings.size();
            CoordMatrix fittingCoords2(4, _numberPoints), otherCoords2(4, _numberPoints);
            int _pointNumber = 0;
            for (auto it = closeMappings.begin(); it != closeMappings.end(); ++it) {
                auto featureNumber = *it;
                auto fittingFeature = fittingFeatures[featureNumber];
                auto otherFeature = featureMap[fittingFeature];
                auto fittingPoint = fittingFeature->getFittingPoint(fittingCoordinates);
                auto otherPoint = otherFeature->getFittingPoint(otherCoordinates);
                for (int i = 0; i < 3; i++) {
                    fittingCoords(i, _pointNumber) = fittingPoint[i];
                    otherCoords(i, _pointNumber) = otherPoint[i];
                }
                fittingCoords(3, _pointNumber) = 1.0;
                otherCoords(3, _pointNumber) = 1.0;
                ++_pointNumber;
            }

            matrix = leastSquaresFit(otherCoords, fittingCoords, boost::none);
        }

        otherCoordinates.transformCoordinates(matrix);
        if (remap) {
            for (size_t i = 0; i < fittingFeatures.size(); i++) {
                if (closeMappings.find(i) == closeMappings.end()) {
                    auto _pos = start + i;
                    integerStringChromosome.setValue(_pos, -1);
                }
                REPORT(Reporter::DEBUG) << "V " << integerStringChromosome.getValue(pos);
            }
        }

        if (Gape::Reporter::isReportingAt(Reporter::DEBUG)) {
            int i = 0;
            for (const auto &mappingFeature: fittingFeatures) {
                auto _pos = start + i;
                REPORT(Reporter::DEBUG) << "V " << integerStringChromosome.getValue(_pos);
                if (auto otherFeature = featureMap.find(mappingFeature); otherFeature != featureMap.end()) {
                    REPORT(Reporter::DEBUG)
                        << mappingFeature->mappingInfo(*otherFeature->second, fittingCoordinates, otherCoordinates);
                }
                i++;
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

        assert(ok && fitted);

        fitness = .0;
        const auto &settings = superpositionGa.getSuperposition().settings.getGapeParameters();
        // Conformational Energies
        if (settings.conformationalWeight > .0) {
            calculateConformationalEnergy();
        }

        // Volume Overlay
        if (settings.volumeWeight > 0) {
            calculateVolumeIntegral();
        }

        // TODO: UserFeatures
        // TODO: Constraints

        // features
        featureOverlay = std::make_shared<FeatureOverlay>(*this);
        const auto featureScore = featureOverlay->scoreOverlay();

        fitness = featureScore
                  + settings.volumeWeight * volumeIntegral
                  - settings.conformationalWeight * conformationalEnergy;

        REPORT(Reporter::DEBUG) << "Fitness " << fitness;

        hasFitness = true;
        return fitness;
        return getFitness();
    }

    bool SuperpositionChromosome::equals(const Gape::SuperpositionChromosome &other) const {
        if (!binaryStringChromosome.equals(other.binaryStringChromosome)) {
            return false;
        }
        return integerStringChromosome.equals(other.integerStringChromosome);
    }

    double SuperpositionChromosome::getFitness() const {
        assert(hasFitness);
        return fitness;
    }

    double SuperpositionChromosome::rebuild() {
        hasFitness = false;
        fitMolecules(false);
        score();
        return getFitness();
    }

    double SuperpositionChromosome::calculateConformationalEnergy() {
        conformationalEnergy = .0;
        const auto &molecules = superpositionGa.getSuperposition().getMolecules();
        conformationalEnergies.clear();
        conformationalEnergies.reserve(molecules.size());
        for (size_t i = 0; i < molecules.size(); i++) {
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
        const auto &molecules = superpositionGa.getSuperposition().getMolecules();
        auto numberMolecules = molecules.size();
        volumeIntegrals = std::make_shared<Array2D<double> >(numberMolecules, numberMolecules);
        volumeIntegral = 0.0;
        double count = 0.0;

        for (size_t i = 0; i < numberMolecules; i++) {
            const auto &moleculeA = molecules[i];
            const auto &conformerA = conformerCoordinates[i]->getConformer();
            for (int j = i + 1; j < numberMolecules; j++) {
                const auto &moleculeB = molecules[j];
                const auto &conformerB = conformerCoordinates[j]->getConformer();
                // we get a cleaner overlay if we compare against everything- at
                // a CPU cost.
                double vol = moleculeA->gaussianIntegral(*moleculeB, conformerA, conformerB);
                vol = (moleculeA->getWeight() + moleculeB->getWeight()) * vol / 2.0;
                volumeIntegrals->set(i, j, vol);
                volumeIntegral += vol;
                REPORT(Reporter::DEBUG) << "Volume Integral for molecules " << i << " and " << j << " " << vol;
                count++;
            }
        }

        volumeIntegral = volumeIntegral / count;
        return volumeIntegral * 2;
    }

    /*
     *
     * Uses distance between features to determine if two chromosomes share a niche.
     *
     */
    bool SuperpositionChromosome::sameNiche(const SuperpositionChromosome &other) const {
        const auto &settings = superpositionGa.getSuperposition().settings.getGapeParameters();
        if (!settings.useNiches) {
            return false;
        }

        const auto featurePointSets = featureOverlay->getFeaturePointSets();
        const auto otherFeaturePointSets = other.featureOverlay->getFeaturePointSets();
        const int numberFeatures = std::accumulate(featurePointSets.begin(), featurePointSets.end(), 0,
                                                   [](const auto &sum, const auto &pair) {
                                                       return sum + pair.second->getFeatures().size();
                                                   });
        const auto nicheDistance = settings.nicheDistance;
        const double sqrSumTolerance = numberFeatures * nicheDistance * nicheDistance;
        double sumTolerance = 0.0;
        for (const auto &[featureType, featurePointSet]: featurePointSets) {
            const auto &features = featurePointSet->getFeatures();
            const auto &otherFeatures = otherFeaturePointSets.at(featureType)->getFeatures();
            assert(features.size() == otherFeatures.size());
            for (int i = 0; i < features.size(); i++) {
                const auto &featurePoint = features.at(i)->point;
                const auto &otherFeaturePoint = otherFeatures.at(i)->point;
                sumTolerance += (featurePoint - otherFeaturePoint).lengthSq();
                if (sumTolerance > sqrSumTolerance) {
                    return true;
                }
            }
        }

        return false;
    }

    std::string SuperpositionChromosome::info() const {
        const auto format = boost::format("Fit %7.1f [Donor %2.1f Acc %2.1f Ring %2.1f Conf %5.1f Vol %5.1f") % fitness
                            % featureOverlay->getDonorHydrogenScore() % featureOverlay->getAcceptorAtomScore() %
                            featureOverlay->getAromaticRingScore() % conformationalEnergy % volumeIntegral;
        return format.str();
    }

    void SuperpositionChromosome::outputSolution(std::ostream &outStream, const std::string &prefix) const {
        SDWriter writer(&outStream);
        std::vector<std::string> props {"CONFORMATIONAL_ENERGY"};
        writer.setProps(props);
        const auto &superposition = superpositionGa.getSuperposition();
        for (int i = 0; i < superpositionGa.numberMolecules(); i++) {
            const auto &mol = superpositionGa.getSuperposition().getMolecules()[i];
            RWMol rdMol(mol->getMol());

            rdMol.clearConformers();
            auto &conformer = fittedCoordinates[i]->getConformer();
            rdMol.addConformer(&conformer);

            rdMol.setProp("CONFORMATIONAL_ENERGY", conformationalEnergies[i]);
            auto name = rdMol.getProp<std::string>(common_properties::_Name);
            if (mol.get() == superposition.getBaseMolecule()) {
                name = prefix + name + " Base Molecule";
            } else if (mol.get() == superposition.getFittingMolecule()) {
                name = prefix + name + " Fitting Molecule";
            } else {
                name = prefix + name;
            }
            rdMol.setProp(common_properties::_Name, name);
            rdMol.setProp("CONFORMATIONAL_ENERGY", conformationalEnergies[i]);
            writer.write(rdMol);
        }
        writer.close();

    }

    void SuperpositionChromosome::outputPharmacophoreAsMol(std::ostream &outStream, const std::string &prefix) const {
        SDWriter writer(&outStream);
        RWMol pharmMol;
        Conformer pharmConformer(featureOverlay->numberPharmacophorePoints());
        pharmMol.setProp("FITNESS", fitness);
        pharmMol.setProp("VOLUME_INTEGRAL", volumeIntegral);
        pharmMol.setProp("CONFORMATIONAL_ENERGY", conformationalEnergy);
        pharmMol.setProp("DONOR_HYDROGEN_SCORE", featureOverlay->getDonorHydrogenScore());
        pharmMol.setProp("ACCEPTOR_ATOM_SCORE", featureOverlay->getAcceptorAtomScore());
        pharmMol.setProp("AROMATIC_RING_SCORE", featureOverlay->getAromaticRingScore());
        int count = 0;
        for (const auto& [featureType, featurePointSet]: featureOverlay->getFeaturePointSets()) {
            for (const auto featurePoint: featurePointSet->getFeatures()) {
                if (featurePoint->isPharmacophoreFeature) {
                    pharmMol.addAtom(new Atom(0));
                    pharmConformer.setAtomPos(count, featurePoint->point);
                    ++count;
                }
            }
        }

        writer.write(pharmMol);
        writer.close();
    }
}
