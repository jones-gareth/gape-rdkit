#include "SuperpositionGa.h"
#include "util/Reporter.h"
#include "SuperpositionChromosome.h"

namespace Gape {

    SuperpositionGa::SuperpositionGa(const Superposition& superposition) : superposition(superposition),
                                                                           integerStringChromosomePolicy(
                                                                               getRng(),
                                                                               superposition.getIntegerStringLength()),
                                                                           binaryStringChromosomePolicy(getRng()) {
        integerStringChromosomePolicy.setAllowNulls(true);
        int pos = 0;
        for (const auto range: superposition.getIntegerStringRanges()) {
            integerStringChromosomePolicy.setMax(pos, range);
            ++pos;
        }

        const auto& gapeParameters = superposition.settings.getGapeParameters();
        int popsize = gapeParameters.populationSize;
        int numberIslands = gapeParameters.numberIslands;
        int numberIterations = gapeParameters.numberIterations;

        if (gapeParameters.guessGaParameters) {
            REPORT(Reporter::INFO) << "Guessing GA Paramteters";
            const auto numberMolecules = static_cast<int>(superposition.getMolecules().size());
            numberIslands = numberMolecules + 1;
            numberIterations = numberMolecules * 15000;
            popsize = 100;
            REPORT(Reporter::INFO) << "Number islands " <<numberIslands << " popsize " << popsize + " number iterations" << numberIterations;
        }

        setPopsize(popsize);
        setSelectionPressure(gapeParameters.selectionPressure);

    }

    std::shared_ptr<SuperpositionChromosome> SuperpositionGa::run(int runNumber) {
        SuperpositionGaPopulation population(*this);
        return nullptr;
    }


    const std::vector<std::shared_ptr<GaOperation<SuperpositionChromosome>>> SuperpositionGa::getOperations() const {

    }

    void SuperpositionGa::superpositionMutateOperation(const std::vector<std::shared_ptr<SuperpositionChromosome>>& parents, std::vector<std::shared_ptr<SuperpositionChromosome>>& children) {

    }

    void SuperpositionGa::superpositionCrossoverOperation(const std::vector<std::shared_ptr<SuperpositionChromosome>>& parents, std::vector<std::shared_ptr<SuperpositionChromosome>>& children) {

    }



}
