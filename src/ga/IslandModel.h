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

    public:
        explicit IslandModel(PopulationPolicy &populationPolicy);

        ~IslandModel() = default;

        void create();

        void iterate();

        void rebuild();

        [[nodiscard]] std::string populationInfo() const;

        [[nodiscard]] std::string info() const;

        [[nodiscard]] double getBestScore() const { return bestScore; }

        std::shared_ptr<Chromosome> getBest();
    };

    template<typename Chromosome, typename PopulationPolicy>
    IslandModel<Chromosome,
        PopulationPolicy>::IslandModel(PopulationPolicy &populationPolicy): populationPolicy(populationPolicy),
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
            if (operation->operationName == OperationName::Migrate) {
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
            if (testFitness > bestScore)
                bestScore = testFitness;
        }
    }

    template<typename Chromosome, typename PopulationPolicy>
    void IslandModel<Chromosome, PopulationPolicy>::iterate() {
        double val = rng.normalRand() * totalOperatorWeights;
        auto &pop = populations[currentPopulationNumber];
        if (val < migrationWeight) {
            auto otherPopulationNumber = rng.randomInt(0, numberIslands-1);
            if (otherPopulationNumber == currentPopulationNumber)
                otherPopulationNumber++;
            auto& parent = populations[otherPopulationNumber]->selectParent();
            pop->addToPopulation(parent);
            numberMigrations++;
        } else {
            pop->iterate();
        }

        numberOperations++;
        auto testFitness = pop->getBestScore();
        if (testFitness > bestScore) {
            testFitness = bestScore;
            const auto format = boost::format("Island Pop %2d Op %5d new best: ") % (currentPopulationNumber+1) % numberOperations;
            REPORT(Reporter::DETAIL) << format << pop->getBest()->info();
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

        for (size_t i = 0; i < numberIslands; i++) {
            ss << "Island " << i + 1 << " " << populations.at(i)->info();
        }
        return ss.str();
    }

    template<typename Chromosome, typename PopulationPolicy>
    shared_ptr<Chromosome> IslandModel<Chromosome, PopulationPolicy>::getBest() {
        auto best = std::max_element(populations.begin(), populations.end(), [](const auto &a, const auto &b) {
            return a->getBestScore() < b->getBestScore();
        });
        const auto c = *best;
        return c->getBest();
    }
}
