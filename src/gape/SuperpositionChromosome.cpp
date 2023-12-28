//
// Created by gareth on 10/18/22.
//

#include "SuperpositionChromosome.h"

#include "SuperpositionGa.h"

namespace Gape {
    SuperpositionChromosome::SuperpositionChromosome(const SuperpositionGa& superpositionGa):
        superpositionGa(superpositionGa),
    integerStringChromosome(superpositionGa.getSuperposition().getIntegerStringLength(), superpositionGa.getRng(), superpositionGa.getIntgerStringChromosomePolicy()),
    binaryStringChromosome(superpositionGa.getSuperposition().getBinaryStringLength(), superpositionGa.getRng(), superpositionGa.getBinaryStringChromosomePolicy())
    {

    }

    void SuperpositionChromosome::copyGene(const Chromosome& c) {
        const auto& parent = dynamic_cast<const SuperpositionChromosome&>(c);
        integerStringChromosome.copyGene(parent.integerStringChromosome);
        binaryStringChromosome.copyGene(parent.binaryStringChromosome);
    }

}
