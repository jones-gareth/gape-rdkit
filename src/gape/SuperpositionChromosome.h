//
// Created by gareth on 10/18/22.
//

#pragma once

#include "SuperpositionGa.h"
#include "ga/StringChromosome.h"

namespace Gape {
    class SuperpositionGa;

    class SuperpositionChromosome {
        OperationName operationName = OperationName::None;
        bool fitted = false;
        std::vector<std::shared_ptr<SuperpositionCoordinates>> conformerCoordinates;
        std::vector<std::shared_ptr<SuperpositionCoordinates>> fittedCoordinates;

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

        double score();

        void fitMolecules(bool remap);

        bool fitMolecule(int start, const SuperpositionMolecule& fittingMolecule,
                         const SuperpositionMolecule& otherMolecule,
                         const SuperpositionCoordinates& fittingCoordinates, SuperpositionCoordinates otherCoordinatesi,
                         bool remap = true);
    };
} // GapeApp
