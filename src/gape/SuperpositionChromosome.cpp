//
// Created by gareth on 10/18/22.
//

#include "SuperpositionChromosome.h"

#include "SuperpositionGa.h"

namespace Gape {
    SuperpositionChromosome::SuperpositionChromosome(const SuperpositionGa& superpositionGa):
        superpositionGa(superpositionGa),
    integerStringChromosome(superpositionGa.getSuperposition().getIntegerStringLength(), superpositionGa.getRng(), superpositionGa.getIntgerStringChromosomePolicy())
    {

    }

}
