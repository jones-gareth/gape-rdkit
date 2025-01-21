//
// Created by gareth on 10/18/22.
//

#pragma once

#include <util/Array2D.h>

#include "SuperpositionGa.h"
#include "ga/StringChromosome.h"

namespace Gape {
    class SuperpositionGa;
    class FeatureOverlay;

    struct FeatureInformation {
        bool isPharmacophoreFeature = false;
        int numberMatches = 0;
        const Feature *feature;
        const RDGeom::Point3D point;
        const SuperpositionCoordinates &coordinates;

        FeatureInformation(const Feature *feature, const SuperpositionCoordinates &coordinates): feature(feature),
            point(feature->getFittingPoint(coordinates)), coordinates(coordinates) {
        }
    };

    class SuperpositionChromosome {
        static const double passOneDistance;
        static const double maxPassDistance;

        OperationName operationName = OperationName::None;
        bool fitted = false;
        bool hasFitness = false;
        bool ok = false;
        std::vector<std::shared_ptr<SuperpositionCoordinates> > conformerCoordinates;
        std::vector<std::shared_ptr<SuperpositionCoordinates> > fittedCoordinates;
        std::vector<double> conformationalEnergies;
        std::shared_ptr<Array2D<double>> volumeIntegrals = nullptr;
        double conformationalEnergy = 0.0;
        double volumeIntegral = 0.0;
        double fitness = 0.0;
        std::shared_ptr<FeatureOverlay> featureOverlay = nullptr;

        double calculateConformationalEnergy();

        double calculateVolumeIntegral();

    public:
        const SuperpositionGa &superpositionGa;
        BinaryStringChromosome binaryStringChromosome;
        IntegerStringChromosome integerStringChromosome;

        SuperpositionChromosome() = delete;

        SuperpositionChromosome(const SuperpositionChromosome &) = delete;

        SuperpositionChromosome &operator=(const SuperpositionChromosome &) = delete;

        explicit SuperpositionChromosome(const SuperpositionGa &superpositionGa);

        void setOperationName(const OperationName operationName_) { operationName = operationName_; }

        const OperationName &getOperationName() const { return operationName; }

        void copyGene(const SuperpositionChromosome &from);

        void initialize();

        double getFitness() const;

        double score();

        void fitMolecules(bool remap = true);

        bool fitMolecule(int start, const SuperpositionMolecule &fittingMolecule,
                         const SuperpositionMolecule &otherMolecule,
                         const SuperpositionCoordinates &fittingCoordinates, SuperpositionCoordinates &otherCoordinates,
                         bool remap = true);

        bool isOk();

        double rebuild();

        bool equals(const SuperpositionChromosome &other) const;

        const std::vector<std::shared_ptr<SuperpositionCoordinates>>& getFittedCoordinates() const {return fittedCoordinates;}

        bool sameNiche(const SuperpositionChromosome &other) const;

        std::string info() const;

        void outputSolution(std::ostream& outStream, const std::string& prefix) const;

        void outputPharmacophoreAsMol(std::ostream& outStream, const std::string &prefix) const;

    };
} // GapeApp
