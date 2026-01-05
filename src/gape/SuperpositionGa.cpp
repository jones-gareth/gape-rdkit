#include "SuperpositionGa.h"

#include <chrono>
#include <future>
#include <GraphMol/FileParsers/MolSupplier.h>

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
            REPORT(Reporter::INFO) << "Number islands " << numberIslands << " popsize " << popsize <<
 " number iterations" << numberOperations;
        }

        setPopsize(popsize);
        setSelectionPressure(gapeParameters.selectionPressure);
    }

    std::shared_ptr<SuperpositionChromosome> SuperpositionGa::run(int runNumber) {
        const auto runStart = std::chrono::high_resolution_clock::now();

        const auto &gapeParameters = superposition.settings.getGapeParameters();

        double startFittingRadius = gapeParameters.startFittingRadius;
        double finishFittingRadius = gapeParameters.finishFittingRadius;
        bool scaleFitting = gapeParameters.scaleFitting;
        if (scaleFitting)
            Feature::setRadius(startFittingRadius);
        else
            Feature::setRadius(finishFittingRadius);

        SuperpositionGaPopulation population(runNumber, *this);
        auto format =
                boost::format(
                    "Running GA run %2d number operations %5d number islands %5d population size %5d ") %
                runNumber % numberOperations % getNumberIslands() % getPopsize();
        REPORT(Reporter::INFO) << format.str();
        population.create();
        REPORT(Reporter::DEBUG) << population.populationInfo();
        double bestScore = population.getBestScore();

        nichesOn = gapeParameters.useNiches;
        double nichingOff = gapeParameters.nichingOff;
        int numberRebuilds = gapeParameters.numberRebuilds;
        int reportInterval = 5000;
        double d = numberOperations * nichingOff;
        int opOff = static_cast<int>(std::ceil(d));
        REPORT(Reporter::INFO) << "Run " << runNumber << " Turning off niching and fitting radius scaling at " << opOff;
        int rebuildInterval = opOff / numberRebuilds;

        REPORT(Reporter::INFO) << population.info() << endl;

        for (int i = 0; i < numberOperations; i++) {
            // Turn off niching after 80% to allow convergence
            if (nichesOn && i >= opOff) {
                REPORT(Reporter::INFO) << "Run " << runNumber << " OP " << i << ": turning off niching";
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
                REPORT(Reporter::DEBUG) << "Run " << runNumber << " setting fitting radius to " << Feature::getRadius();
                population.rebuild();
            }

            population.iterate();

            if (reportInterval > 0 && i > 1 && (i + 1) % reportInterval == 0)
                REPORT(Reporter::INFO) << population.info();
        }

        const auto best = population.getBest();
        best->solutionNumber = runNumber;

        const auto n = runNumber + 1;
        REPORT(Reporter::INFO) << "Best solution for run number " << n << ":" << best->info();
        const auto fileName = boost::format("GA_solution_%d.sdf") % n;
        const auto prefix = boost::format("Solution number %d") % n;
        std::fstream out{fileName.str(), std::fstream::out};
        best->outputSolution(out, prefix.str());
        out.close();
        const auto runEnd = std::chrono::high_resolution_clock::now();
        const auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(runEnd - runStart);
        const auto timeString = boost::format("Run time %5.2f") % (duration.count() / 1000.0);
        REPORT(Reporter::INFO) << timeString.str();
        return best;
    }

    std::vector<std::shared_ptr<GaOperation<SuperpositionChromosome> > > SuperpositionGa::getOperations() const {
        const auto &parameters = superposition.settings.getGapeParameters();
        const auto mutationOperation = std::make_shared<GaOperation<SuperpositionChromosome> >(OperationName::Mutate,
            1, 1, parameters.mutationWeight, &superpositionMutateOperation);
        const auto crossoverOperation = std::make_shared<GaOperation<SuperpositionChromosome> >(
            OperationName::Crossover,
            2, 2, parameters.crossoverWeight, &superpositionCrossoverOperation);
        const auto migrationOperation = std::make_shared<GaOperation<SuperpositionChromosome> >(
            OperationName::Migrate, 1, 1, parameters.migrationWeight, &superpositionMigrationOperation);
        std::vector operations{mutationOperation, crossoverOperation, migrationOperation};
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
            child->copyConformerCoordinates(*parent);
            child->integerStringChromosome.mutate();
        }
    }

    void SuperpositionGa::superpositionCrossoverOperation(
        const std::vector<std::shared_ptr<SuperpositionChromosome> > &parents,
        std::vector<std::shared_ptr<SuperpositionChromosome> > &children) {
        assert(parents.size() == 2);
        assert(children.size() == 2);
        children[0]->copyGene(*parents[0]);
        children[1]->copyGene(*parents[1]);
        if (const auto &parent1BinaryString = parents[0]->binaryStringChromosome;
            parent1BinaryString.getLength() > 0 && parent1BinaryString.getRng().randomBoolean()) {
            children[0]->setOperationName(OperationName::BinaryStringCrossover);
            children[1]->setOperationName(OperationName::BinaryStringCrossover);
            parent1BinaryString.onePointCrossover(parents[1]->binaryStringChromosome,
                                                  children[0]->binaryStringChromosome,
                                                  children[1]->binaryStringChromosome);
        } else {
            children[0]->copyConformerCoordinates(*parents[0]);
            children[1]->copyConformerCoordinates(*parents[1]);
            children[0]->setOperationName(OperationName::IntegerStringCrossover);
            children[1]->setOperationName(OperationName::IntegerStringCrossover);
            parents[0]->integerStringChromosome.fullMixing(parents[1]->integerStringChromosome,
                                                           children[0]->integerStringChromosome,
                                                           children[1]->integerStringChromosome);
        }
    }

    void SuperpositionGa::superpositionMigrationOperation(
        // this function is never called migration is handled in IslandModel::iterate
        const std::vector<std::shared_ptr<SuperpositionChromosome> > &parents,
        std::vector<std::shared_ptr<SuperpositionChromosome> > &children) {
        assert(parents.size() == 1);
        assert(children.size() == 1);
        children[0]->copyGene(*parents[0]);
        children[0]->setOperationName(OperationName::Migrate);
    }

    std::shared_ptr<SuperpositionChromosome> SuperpositionGa::singleRun(const Superposition &superposition,
                                                                        const int runNumber) {
        SuperpositionGa ga(superposition);
        return ga.run(runNumber);
    }

    std::vector<std::shared_ptr<SuperpositionChromosome> >
    SuperpositionGa::batchRun(const Superposition &superposition) {
        std::vector<std::shared_ptr<SuperpositionChromosome> > solutions;
        const auto batchStart = std::chrono::high_resolution_clock::now();
        const auto numberRuns = superposition.settings.getGapeParameters().numberRuns;
        solutions.clear();
        solutions.reserve(numberRuns);
        // bool parallelRuns = superposition.settings.getGapeParameters().parallelRuns;
        bool parallelRuns = true;
        if (parallelRuns) {
            /*
            std::vector<future<std::shared_ptr<SuperpositionChromosome>>> tasks;
            tasks.reserve(numberRuns);
            for (int runNumber = 0; runNumber < numberRuns; runNumber++) {
                auto future = std::async(std::launch::async, &Gape::SuperpositionGa::singleRun, std::ref(superposition), runNumber);
                tasks.push_back(std::move(future));
            }
            std::transform(tasks.begin(), tasks.end(), back_inserter(solutions),
               [](future<std::shared_ptr<SuperpositionChromosome>> &f) { return f.get(); });
               */
            solutions.resize(numberRuns);
            vector<std::unique_ptr<std::thread> > threads;
            threads.reserve(numberRuns);
            for (int runNumber = 0; runNumber < numberRuns; runNumber++) {
                auto lambda = [&superposition, runNumber, &solutions]() {
                    const auto best = singleRun(superposition, runNumber);
                    solutions[runNumber] = best;
                };
                auto t = std::make_unique<std::thread>(lambda);
                threads.push_back(std::move(t));
            }
            for (auto &t: threads) {
                t->join();
            }
        } else {
            for (int runNumber = 0; runNumber < numberRuns; runNumber++) {
                const auto best = singleRun(superposition, runNumber);
                solutions.push_back(best);
            }
        }

        const auto batchEnd = std::chrono::high_resolution_clock::now();
        const auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(batchEnd - batchStart);
        const auto timeString = boost::format("Batch time %5.2f") % (duration.count() / 1000.0);
        REPORT(Reporter::INFO) << timeString.str();

        if (numberRuns > 0) {
            std::sort(solutions.begin(), solutions.end(),
                      [](const auto &a, const auto &b) {
                          return a->getFitness() > b->getFitness();
                      });


            for (int runNumber = 0; runNumber < numberRuns; runNumber++) {
                const auto &c = solutions[runNumber];
                const auto n = runNumber + 1;
                const auto msg = boost::format("Rank %2d solution %2d ") % n % (c->solutionNumber + 1);
                REPORT(Reporter::INFO) << msg.str() << c->info();
                const auto fileName = boost::format("GA_rank_%d.sdf") % n;
                const auto prefix = boost::format("Solution rank %d") % n;
                std::fstream out{fileName.str(), std::fstream::out};
                c->outputSolution(out, prefix.str());
                out.close();
            }
        }

        return solutions;
    }

    std::vector<std::shared_ptr<SuperpositionChromosome> >
    SuperpositionGa::testBatchRun(const std::string &inputFile, const GapeSettings &settings) {
        std::vector<std::shared_ptr<SuperpositionChromosome> > solutions;
        const auto batchStart = std::chrono::high_resolution_clock::now();
        const auto numberRuns = settings.getGapeParameters().numberRuns;
        solutions.clear();
        solutions.reserve(numberRuns);

        solutions.resize(numberRuns);
        vector<std::unique_ptr<std::thread> > threads;
        threads.reserve(numberRuns);
        for (int runNumber = 0; runNumber < numberRuns; runNumber++) {
            auto lambda = [&inputFile, &settings, &solutions, runNumber]() {
                RDKit::SmilesMolSupplier smilesMolSupplier(inputFile, " ", 0, 1, false, true);
                auto molecules = SuperpositionMolecule::loadMolecules(smilesMolSupplier, settings);
                Superposition superposition(molecules, settings);
                const auto best = singleRun(superposition, runNumber);
                solutions[runNumber] = best;
            };
            auto t = std::make_unique<std::thread>(lambda);
            threads.push_back(std::move(t));
        }
        for (auto &t: threads) {
            t->join();
        }

        const auto batchEnd = std::chrono::high_resolution_clock::now();
        const auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(batchEnd - batchStart);
        const auto timeString = boost::format("Batch time %5.2f") % (duration.count() / 1000.0);
        REPORT(Reporter::INFO) << timeString.str();

        if (numberRuns > 0) {
            std::sort(solutions.begin(), solutions.end(),
                      [](const auto &a, const auto &b) {
                          return a->getFitness() > b->getFitness();
                      });


            for (int runNumber = 0; runNumber < numberRuns; runNumber++) {
                const auto &c = solutions[runNumber];
                const auto n = runNumber + 1;
                const auto msg = boost::format("Rank %2d solution %2d ") % n % (c->solutionNumber + 1);
                REPORT(Reporter::INFO) << msg.str() << c->info();
                const auto fileName = boost::format("GA_rank_%d.sdf") % n;
                const auto prefix = boost::format("Solution rank %d") % n;
                std::fstream out{fileName.str(), std::fstream::out};
                c->outputSolution(out, prefix.str());
                out.close();
            }
        }

        return solutions;
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
