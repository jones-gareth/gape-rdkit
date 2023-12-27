#pragma once

#include "Superposition.h"

#include <ga/GaBase.h>
#include <ga/LinkedPopLinearSel.h>

#include "ga/BinaryStringChromosomePolicy.h"
#include "ga/IntegerStringChromosomePolicy.h"

namespace Gape {
    class SuperpositionChromosome;
    class SuperpositionGa;

    typedef LinkedPopLinearSel<SuperpositionChromosome, SuperpositionGa>
    SuperpositionGaPopulation;

    typedef enum {
        RgroupMutate = 0x01,
        Crossover = 0x02,
        Create = 0x04,
        Migrate = 0x08,
    } OperationName;

    class SuperpositionGa : public GaBase {
        const Superposition& superposition;
        IntegerStringChromosomePolicy integerStringChromosomePolicy;
        BinaryStringChromosomePolicy binaryStringChromosomePolicy;

    public:
        explicit SuperpositionGa(const Superposition& superposition);

        SuperpositionGa(const SuperpositionGa&) = delete;

        SuperpositionGa& operator=(const SuperpositionGa&) = delete;

        std::shared_ptr<SuperpositionChromosome> run(int runNumber);

        [[nodiscard]] const std::vector<std::shared_ptr<GaOperation<SuperpositionChromosome>>>
        getOperations() const;

        static void superpositionMutateOperation(
            const std::vector<std::shared_ptr<SuperpositionChromosome>>&
            parents,
            std::vector<std::shared_ptr<SuperpositionChromosome>>& children);

        static void superpositionCrossoverOperation(
            const std::vector<std::shared_ptr<SuperpositionChromosome>>&
            parents,
            std::vector<std::shared_ptr<SuperpositionChromosome>>& children);

        [[nodiscard]] const Superposition& getSuperposition() const { return superposition; }

        [[nodiscard]] const IntegerStringChromosomePolicy& getIntgerStringChromosomePolicy() const {
            return integerStringChromosomePolicy;
        }

        [[nodiscard]] const BinaryStringChromosomePolicy& getBinaryStringChromosomePolicy() const {
            return binaryStringChromosomePolicy;
        }
    };
}
