#pragma once

#include "Superposition.h"

#include <ga/GaBase.h>
#include <ga/LinkedPopLinearSel.h>

#include "ga/BinaryStringChromosomePolicy.h"
#include "ga/IntegerStringChromosomePolicy.h"
#include "ga/IslandModel.h"

namespace Gape {
    class SuperpositionChromosome;
    class SuperpositionGa;

    typedef IslandModel<SuperpositionChromosome, SuperpositionGa>
    SuperpositionGaPopulation;


    class SuperpositionGa : public GaBase {
        const Superposition& superposition;
        IntegerStringChromosomePolicy integerStringChromosomePolicy;
        BinaryStringChromosomePolicy binaryStringChromosomePolicy;
        int numberOperations = 0;
        bool nichesOn = false;
        std::vector<std::shared_ptr<SuperpositionChromosome>> solutions;

        std::shared_ptr<SuperpositionChromosome> run(int runNumber);

    public:
        explicit SuperpositionGa(const Superposition& superposition);

        SuperpositionGa(const SuperpositionGa&) = delete;

        SuperpositionGa& operator=(const SuperpositionGa&) = delete;

        [[nodiscard]] std::vector<std::shared_ptr<GaOperation<SuperpositionChromosome>>>
        getOperations() const;

        static void superpositionMutateOperation(
            const std::vector<std::shared_ptr<SuperpositionChromosome>>&
            parents,
            std::vector<std::shared_ptr<SuperpositionChromosome>>& children);

        static void superpositionCrossoverOperation(
            const std::vector<std::shared_ptr<SuperpositionChromosome>>&
            parents,
            std::vector<std::shared_ptr<SuperpositionChromosome>>& children);

        static void superpositionMigrationOperation(
               const std::vector<std::shared_ptr<SuperpositionChromosome>>&
               parents,
               std::vector<std::shared_ptr<SuperpositionChromosome>>& children);

        [[nodiscard]] const Superposition& getSuperposition() const { return superposition; }

        [[nodiscard]] const IntegerStringChromosomePolicy& getIntegerStringChromosomePolicy() const {
            return integerStringChromosomePolicy;
        }

        [[nodiscard]] const BinaryStringChromosomePolicy& getBinaryStringChromosomePolicy() const {
            return binaryStringChromosomePolicy;
        }

        int numberMolecules() const {
            return static_cast<int>(superposition.getMolecules().size());
        }

        void run();

        std::shared_ptr<SuperpositionChromosome> createChromosome();

        [[nodiscard]] bool useNiches() const;

        [[nodiscard]] int getNicheSize() const;

        [[nodiscard]] int getNumberIslands() const;
    };
}
