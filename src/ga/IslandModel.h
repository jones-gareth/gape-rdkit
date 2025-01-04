//
// Created by jones on 1/3/2025.
//

#pragma once
#include <memory>
#include <util/RandomUtil.h>

#include "GaOperation.h"
#include "LinkedPopLinearSel.h"

namespace gape {
    template<typename Chromosome, typename PopulationPolicy>
    class IslandModel {
        PopulationPolicy &populationPolicy;
        RandomUtil &rng;
        const size_t numberIslands;
        const std::vector<std::shared_ptr<Gape::GaOperation<Chromosome> > > operations;

        std::vector<Gape::LinkedPopLinearSel<Chromosome, PopulationPolicy> > populations;

    public:
        IslandModel(PopulationPolicy populationPolicy);
    };

    template<typename Chromosome, typename PopulationPolicy>
    IslandModel<Chromosome,
        PopulationPolicy>::IslandModel(PopulationPolicy populationPolicy): populationPolicy(populationPolicy_),
                                                                           rng(populationPolicy.getRng()),
                                                                           numberIslands(
                                                                               populationPolicy.getNumberIslands()) {
        for (size_t i = 0; i < numberIslands; i++) {
            populations.emplace_back(populationPolicy);
        }

    }
}
