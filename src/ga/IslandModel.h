//
// Created by jones on 1/3/2025.
//

#pragma once
#include <memory>
#include <util/RandomUtil.h>

#include "GaOperation.h"
#include "LinkedPopLinearSel.h"

namespace Gape {
    template<typename Chromosome, typename PopulationPolicy>
    class IslandModel {
        const size_t runNumber;
        PopulationPolicy &populationPolicy;
        RandomUtil &rng;
        const size_t numberIslands;
        std::vector<std::shared_ptr<LinkedPopLinearSel<Chromosome, PopulationPolicy> > > populations;
        int currentPopulationNumber = 0;
        double totalOperatorWeights = 0;
        double migrationWeight = 0;
        const std::vector<std::shared_ptr<GaOperation<Chromosome> > > operations;
        std::shared_ptr<GaOperation<Chromosome> > migrationOperation;
        size_t numberOperations = 0;
        size_t numberMigrations = 0;
        double bestScore = -std::numeric_limits<double>::max();
        std::shared_ptr<Chromosome> bestChromosome;

    public:
        explicit IslandModel(size_t runNumber, PopulationPolicy &populationPolicy);

        ~IslandModel() = default;

        void create();

        void iterate();

        void rebuild();

        [[nodiscard]] std::string populationInfo() const;

        [[nodiscard]] std::string info() const;

        [[nodiscard]] double getBestScore() const { return bestScore; }

        std::shared_ptr<Chromosome> getBestChromosome() const { return bestChromosome; }

        std::shared_ptr<Chromosome> getBest();
    };

    template<typename Chromosome, typename PopulationPolicy>
    IslandModel<Chromosome,
        PopulationPolicy>::IslandModel(const size_t runNumber,
                                       PopulationPolicy &populationPolicy) : runNumber(runNumber),
                                                                             populationPolicy(populationPolicy),
                                                                             rng(populationPolicy.getRng()),
                                                                             numberIslands(
                                                                                 populationPolicy.getNumberIslands()),
                                                                             operations(
                                                                                 populationPolicy.getOperations()) {
        populations.reserve(populationPolicy.getNumberIslands());
        for (size_t i = 0; i < numberIslands; i++) {
            populations.
                    push_back(std::make_shared<LinkedPopLinearSel<Chromosome, PopulationPolicy> >(populationPolicy));
        }

        totalOperatorWeights = 0;
        for (auto &operation: operations) {
            totalOperatorWeights += operation->getWeight();
            if (operation->operationName == OperationName::Migrate && numberIslands > 1) {
                migrationWeight += operation->getWeight();
                migrationOperation = operation;
            }
        }
    }

    template<typename Chromosome, typename PopulationPolicy>
    void IslandModel<Chromosome, PopulationPolicy>::create() {
        bestScore = -std::numeric_limits<double>::max();
        for (auto population: populations) {
            population->create();
            auto testFitness = population->getBestScore();
            if (testFitness > bestScore) {
                bestScore = testFitness;
                bestChromosome = population->getBest();
            }
        }
    }

    template<typename Chromosome, typename PopulationPolicy>
    void IslandModel<Chromosome, PopulationPolicy>::iterate() {
        double val = rng.normalRand() * totalOperatorWeights;
        auto &pop = populations[currentPopulationNumber];
        if (val < migrationWeight) {
            auto otherPopulationNumber = rng.randomInt(0, static_cast<int>(numberIslands) - 1);
            if (otherPopulationNumber == currentPopulationNumber)
                otherPopulationNumber++;
            auto &parent = populations[otherPopulationNumber]->selectParent();
            auto child = pop->fetchChild();

            child->copyGene(*parent);
            child->setOperationName(OperationName::Migrate);
            if (child->isOk()) {
                child->score();
                pop->addToPopulation(child);
                numberMigrations++;
            }
        } else {
            pop->iterate();
        }

        numberOperations++;
        auto testFitness = pop->getBestScore();
        if (testFitness > bestScore) {
            testFitness = bestScore;
            bestChromosome = pop->getBest();
            const auto format = boost::format("Island Pop %2d Op %5d Mig %5d new best: ") % (
                                    currentPopulationNumber + 1) % numberOperations % numberMigrations;
            REPORT(Reporter::DEBUG) << format << bestChromosome->info();
        }
        currentPopulationNumber++;
        if (currentPopulationNumber == numberIslands)
            currentPopulationNumber = 0;
    }

    template<typename Chromosome, typename PopulationPolicy>
    void IslandModel<Chromosome, PopulationPolicy>::rebuild() {
        // TODO: chromosomes created by migrations may be shared between populations and we don't want to
        // score them twice (which we currently do)
        for (auto population: populations) {
            population->rebuild();
        }
    }

    template<typename Chromosome, typename PopulationPolicy>
    std::string IslandModel<Chromosome, PopulationPolicy>::populationInfo() const {
        std::stringstream ss;

        for (size_t i = 0; i < numberIslands; i++) {
            ss << "Island " << i + 1 << std::endl;
            ss << populations.at(i)->populationInfo();
        }
        return ss.str();
    }

    template<typename Chromosome, typename PopulationPolicy>
    std::string IslandModel<Chromosome, PopulationPolicy>::info() const {
        std::stringstream ss;

        const auto format = boost::format("Run %2d, Op %5d Mig %5d \nBest: ") % runNumber % numberOperations %
                            numberMigrations;
        ss << format << bestChromosome->info() << std::endl;
        for (size_t i = 0; i < numberIslands; i++) {
            ss << "Island " << i + 1 << " " << populations.at(i)->info() << std::endl;
        }
        return ss.str();
    }

    template<typename Chromosome, typename PopulationPolicy>
    shared_ptr<Chromosome> IslandModel<Chromosome, PopulationPolicy>::getBest() {
        auto best = std::max_element(populations.cbegin(), populations.cend(), [](const auto &a, const auto &b) {
            return a->getBestScore() < b->getBestScore();
        });
        const auto c = *best;
        return c->getBest();
    }
}
