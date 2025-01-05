#include "SuperpositionGa.h"
#include "util/Reporter.h"
#include "SuperpositionChromosome.h"

namespace Gape {
    SuperpositionGa::SuperpositionGa(const Superposition &superposition) : superposition(superposition),
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
        binaryStringChromosomePolicy.setAllowSwitch(true);

        const auto &gapeParameters = superposition.settings.getGapeParameters();
        int popsize = gapeParameters.populationSize;
        int numberIslands = gapeParameters.numberIslands;
        numberOperations = gapeParameters.numberOperations;

        if (gapeParameters.guessGaParameters) {
            REPORT(Reporter::INFO) << "Guessing GA Parameters";
            const auto numberMolecules = static_cast<int>(superposition.getMolecules().size());
            numberIslands = numberMolecules + 1;
            numberOperations = numberMolecules * 15000;
            popsize = 100;
            REPORT(Reporter::INFO) << "Number islands " << numberIslands << " popsize " << popsize +
 " number iterations" << numberOperations;
        }

        setPopsize(popsize);
        setSelectionPressure(gapeParameters.selectionPressure);
    }

    std::shared_ptr<SuperpositionChromosome> SuperpositionGa::run(int runNumber) {
        const auto &gapeParameters = superposition.settings.getGapeParameters();

        double startFittingRadius = gapeParameters.startFittingRadius;
        double finishFittingRadius = gapeParameters.finishFittingRadius;
        bool scaleFitting = gapeParameters.scaleFitting;
        if (scaleFitting)
            Feature::setRadius(startFittingRadius);
        else
            Feature::setRadius(finishFittingRadius);

        SuperpositionGaPopulation population(*this);
        auto format =
                boost::format(
                    "Running GA run %2d number operations %5d population size %5d ") %
                runNumber % numberOperations % getPopsize();
        REPORT(Reporter::INFO) << format.str();
        population.create();
        REPORT(Reporter::DEBUG) << population.populationInfo();
        double bestScore = population.getBestScore();

        nichesOn = gapeParameters.useNiches;
        double nichingOff = gapeParameters.nichingOff;
        int numberRebuilds = gapeParameters.numberRebuilds;
        int reportInterval = 1000;
        double d = numberOperations * nichingOff;
        int opOff = static_cast<int>(std::ceil(d));
		REPORT(Reporter::INFO) << "Turning off niching and fitting radius scaling at " << opOff;
		int rebuildInterval = opOff / numberRebuilds;

        REPORT(Reporter::INFO) << population.info() << endl;

        for (int i = 0; i < numberOperations; i++) {
            // Turn off niching after 80% to allow convergence
            if (nichesOn && i >= opOff) {
                REPORT(Reporter::INFO) << "OP " << i << ": turning off niching";
                nichesOn = false;
            }

            // anneal fitting point radius
            if (scaleFitting && i > 0 && i % rebuildInterval == 0) {
                if (i >= opOff) {
                    Feature::setRadius(finishFittingRadius);
                    scaleFitting = false;
                } else {
                    const double frac = static_cast<double>(opOff - i) / static_cast<double>(opOff);
                    const double radius = finishFittingRadius
                            + (startFittingRadius - finishFittingRadius) * frac;
                    Feature::setRadius(radius);
                }
                REPORT(Reporter::INFO) << "Setting fitting radius to " << Feature::getRadius();
                population.rebuild();
            }

            population.iterate();

            if (reportInterval > 0 && i > 1 && (i + 1) % reportInterval == 0)
                REPORT(Reporter::INFO) << population.info();

        }

        const auto best = population.getBest();
        return best;
    }

    std::vector<std::shared_ptr<GaOperation<SuperpositionChromosome> > > SuperpositionGa::getOperations() const {
        const auto &parameters = superposition.settings.getGapeParameters();
        const auto mutationOperation = std::make_shared<GaOperation<SuperpositionChromosome> >(OperationName::Mutate,
            1, 1, parameters.mutationWeight, &superpositionMutateOperation);
        const auto crossoverOperation = std::make_shared<GaOperation<SuperpositionChromosome> >(OperationName::Crossover,
            2, 2, parameters.crossoverWeight, &superpositionCrossoverOperation);
        const auto migrationOperation = std::make_shared<GaOperation<SuperpositionChromosome> >(OperationName::Migrate, 1, 1, parameters.migrationWeight, &superpositionMigrationOperation);
        std::vector operations{mutationOperation, crossoverOperation};
        return operations;
    }

    void SuperpositionGa::superpositionMutateOperation(
        const std::vector<std::shared_ptr<SuperpositionChromosome> > &parents,
        std::vector<std::shared_ptr<SuperpositionChromosome> > &children) {
        assert(parents.size() == 1);
        assert(children.size() == 1);
        auto &parent = parents[0];
        auto &child = children[0];
        child->copyGene(*parent);

        if (auto &childBinaryString = child->binaryStringChromosome;
            childBinaryString.getLength() > 0 && childBinaryString.getRng().randomBoolean()) {
            child->setOperationName(OperationName::BinaryStringMutate);
            childBinaryString.mutate();
        } else {
            child->setOperationName(OperationName::IntegerStringMutate);
            child->integerStringChromosome.mutate();
        }
    }

    void SuperpositionGa::superpositionCrossoverOperation(
        const std::vector<std::shared_ptr<SuperpositionChromosome> > &parents,
        std::vector<std::shared_ptr<SuperpositionChromosome> > &children) {
        assert(parents.size() == 2);
        assert(children.size() == 2);
        if (const auto &parent1BinaryString = parents[0]->binaryStringChromosome;
            parent1BinaryString.getLength() > 0 && parent1BinaryString.getRng().randomBoolean()) {
            children[0]->setOperationName(OperationName::BinaryStringCrossover);
            children[1]->setOperationName(OperationName::BinaryStringCrossover);
            parent1BinaryString.onePointCrossover(parents[1]->binaryStringChromosome,
                                                  children[0]->binaryStringChromosome,
                                                  children[1]->binaryStringChromosome);
        } else {
            children[0]->setOperationName(OperationName::IntegerStringCrossover);
            children[1]->setOperationName(OperationName::IntegerStringCrossover);
            parents[0]->integerStringChromosome.fullMixing(parents[1]->integerStringChromosome,
                                                           children[0]->integerStringChromosome,
                                                           children[1]->integerStringChromosome);
        }
    }

    void SuperpositionGa::superpositionMigrationOperation(
       const std::vector<std::shared_ptr<SuperpositionChromosome> > &parents,
       std::vector<std::shared_ptr<SuperpositionChromosome> > &children) {
        assert(parents.size() == 1);
        assert(children.size() == 1);
        // migration does nothing
    }

    void SuperpositionGa::run() {
        const auto numberRuns = superposition.settings.getGapeParameters().numberRuns;
        for (int runNumber = 0; runNumber < numberRuns; runNumber++) {
            run(runNumber);
        }
    }

    std::shared_ptr<SuperpositionChromosome> SuperpositionGa::createChromosome() {
        return std::make_shared<SuperpositionChromosome>(*this);
    }


    bool SuperpositionGa::useNiches() const {
        if (!nichesOn) return false;
        const auto &settings = getSuperposition().settings.getGapeParameters();
        return settings.useNiches;
    }

    int SuperpositionGa::getNicheSize() const {
        const auto &settings = getSuperposition().settings.getGapeParameters();
        return settings.nicheSize;
    }

    int SuperpositionGa::getNumberIslands() const {
        const auto &settings = getSuperposition().settings.getGapeParameters();
        return settings.numberIslands;
    }
}
