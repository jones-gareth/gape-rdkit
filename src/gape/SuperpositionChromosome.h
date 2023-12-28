//
// Created by gareth on 10/18/22.
//

#pragma once

#include "SuperpositionGa.h"
#include "ga/Chromosome.h"
#include "ga/StringChromosome.h"

namespace Gape {
    class SuperpositionGa;
    enum class OperationName;

    class SuperpositionChromosome : public Chromosome {
        OperationName operationName;

    public:
        const SuperpositionGa& superpositionGa;
        BinaryStringChromosome binaryStringChromosome;
        IntegerStringChromosome integerStringChromosome;


        SuperpositionChromosome() = delete;

        SuperpositionChromosome(const SuperpositionChromosome&) = delete;

        SuperpositionChromosome& operator=(const SuperpositionChromosome&) = delete;

        explicit SuperpositionChromosome(const SuperpositionGa &superpositionGa);

        void setOperationName(const OperationName operationName_) {operationName = operationName_;}

        const OperationName& getOperationName() const {return operationName; }

        Chromosome& create() override;

        double rebuild(Chromosome& c) override;

        bool equals(const Chromosome& c) const override;

        double distance(const Chromosome& c) const override;

        bool sameNiche(const Chromosome& c) const override;

        bool ok() const override;

        void copyGene(const Chromosome& c) override;

        std::string fitnessInfo() const override;

        std::string geneInfo() const override;

        void calculateFitness() override;
    };
} // GapeApp
