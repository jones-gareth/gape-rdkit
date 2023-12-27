//
// Created by gareth on 10/18/22.
//

#pragma once

#include "ga/Chromosome.h"
#include "ga/StringChromosome.h"

namespace Gape {
    class SuperpositionGa;

    class SuperpositionChromosome : public Chromosome {
    public:
        const SuperpositionGa& superpositionGa;
        BinaryStringChromosome binaryStringChromosome;
        IntegerStringChromosome integerStringChromosome;


        SuperpositionChromosome() = delete;

        SuperpositionChromosome(const SuperpositionChromosome&) = delete;

        SuperpositionChromosome& operator=(const SuperpositionChromosome&) = delete;

        explicit SuperpositionChromosome(const SuperpositionGa &superpositionGa);
    };
} // GapeApp
