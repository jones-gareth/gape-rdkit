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
        std::vector<std::shared_ptr<LinkedPopLinearSel<Chromosome, PopulationPolicy> >> populations;
        int currentPopulationNumber = 0;

    public:
        explicit IslandModel(PopulationPolicy &populationPolicy);

        ~IslandModel() = default;

        void create();

        void iterate();

        void rebuild();

        [[nodiscard]] std::string populationInfo() const;

        [[nodiscard]] std::string info() const;

        [[nodiscard]] double getBestScore() const;
    };

    template<typename Chromosome, typename PopulationPolicy>
    IslandModel<Chromosome,
        PopulationPolicy>::IslandModel(PopulationPolicy &populationPolicy): populationPolicy(populationPolicy),
                                                                            rng(populationPolicy.getRng()),
                                                                            numberIslands(
                                                                                populationPolicy.getNumberIslands()) {
        populations.reserve(populationPolicy.getNumberIslands());
        for (size_t i = 0; i < numberIslands; i++) {
            populations.push_back(std::make_shared<LinkedPopLinearSel<Chromosome, PopulationPolicy>>(populationPolicy));
        }
    }

    template<typename Chromosome, typename PopulationPolicy>
    void IslandModel<Chromosome, PopulationPolicy>::create() {
        for (auto population : populations) {
            population->create();
        }
    }

    template<typename Chromosome, typename PopulationPolicy>
    void IslandModel<Chromosome, PopulationPolicy>::iterate() {
        auto& pop = populations[currentPopulationNumber];
        /*
        int wt = rng.randomInt(0, totalOperatorWt + migrateWt - 1);
        if (wt < migrateWt && nIslands > 1) {
            migrate();
        } else {
            Operator op = pop.selectOperator();
            pop.applyOperator(op);
        }
        if (pop.getBest().getFitness() > bestFitness) {
            getBest();
            infoMessageLogger.infoMessage(bestInfo());
        }
        nOperations++;
        currentPopNo++;
        if (currentPopNo == nIslands)
            currentPopNo = 0;
    */
    }

    template<typename Chromosome, typename PopulationPolicy>
    void IslandModel<Chromosome, PopulationPolicy>::rebuild() {
        // TODO: chromosomes created by migrations may be shared between populations and we don't want to
        // score them twice (which we currently do)
        for (auto population : populations) {
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
    double IslandModel<Chromosome, PopulationPolicy>::getBestScore() const {
        double bestScore = -std::numeric_limits<double>::max();
        for (const auto &population : populations) {
            double score = population->getBestScore();
            if (score > bestScore) {
                bestScore = score;
            }
        }
        return bestScore;
    }
}
