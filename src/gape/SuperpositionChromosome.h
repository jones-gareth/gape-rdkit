//
// Created by gareth on 10/18/22.
//

#pragma once

#include "SuperpositionGa.h"
#include "ga/StringChromosome.h"

namespace Gape {
    class SuperpositionGa;

    struct FeatureInformation {
        bool isPharmacophoreFeature = false;
        int numberMatches = 0;
        const Feature * feature;
        const RDGeom::Point3D point;

        FeatureInformation(const Feature * feature, const RDGeom::Point3D & point): feature(feature), point(point) {}
    };

    class SuperpositionChromosome {

        static const double passOneDistance;
        static const double maxPassDistance;

        OperationName operationName = OperationName::None;
        bool fitted = false;
        bool hasFitness = false;
        bool ok = false;
        std::vector<std::shared_ptr<SuperpositionCoordinates>> conformerCoordinates;
        std::vector<std::shared_ptr<SuperpositionCoordinates>> fittedCoordinates;
        std::vector<double> conformationalEnergies;

        std::map<Feature *, std::shared_ptr<FeatureInformation>> featureMapping;
        double fitness;

        double calculateConformationalEnergy();

        double calculateVolumeIntegral();

    public:
        const SuperpositionGa& superpositionGa;
        BinaryStringChromosome binaryStringChromosome;
        IntegerStringChromosome integerStringChromosome;

        SuperpositionChromosome() = delete;

        SuperpositionChromosome(const SuperpositionChromosome&) = delete;

        SuperpositionChromosome& operator=(const SuperpositionChromosome&) = delete;

        explicit SuperpositionChromosome(const SuperpositionGa& superpositionGa);

        void setOperationName(const OperationName operationName_) { operationName = operationName_; }

        const OperationName& getOperationName() const { return operationName; }

        void copyGene(const SuperpositionChromosome& from);

        void initialize();

        double getFitness();

        double score();

        void fitMolecules(bool remap = true);

        bool fitMolecule(int start, const SuperpositionMolecule& fittingMolecule,
                         const SuperpositionMolecule& otherMolecule,
                         const SuperpositionCoordinates& fittingCoordinates, SuperpositionCoordinates &otherCoordinates,
                         bool remap = true);

        bool isOk();

        double rebuild();

        bool equals(const SuperpositionChromosome& other) const;

    };
} // GapeApp
