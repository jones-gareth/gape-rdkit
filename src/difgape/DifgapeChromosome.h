/*
 * DifgapeChromosome.h
 *
 *  Created on: Jun 11, 2014
 *      Author: Gareth Jones
 */

#ifndef DIFGAPECHROMOSOME_H_
#define DIFGAPECHROMOSOME_H_

#include <string>
#include "StringChromosome.h"

namespace Difgape {

// use forward declaration- don't include conformer matcher here
class ConformerMatcher;

using namespace std;
using namespace GapeGa;

class DifgapeChromosome: public IntegerStringChromosome {
public:
	DifgapeChromosome(ConformerMatcher & cm_) ;

	virtual ~DifgapeChromosome() {
		;
	}

	double getFitness() const {
		return fitness;
	}

	bool isOk() {
		return true;
	}

	/**
	 * Calculates the fitness of the selected conformers
	 * @return
	 */
	double score();

	/**
	 * Creates an informational string
	 *
	 * @return
	 */
	 const std::string info() const;

	/**
	 * Creates a full summary string
	 *
	 * @return
	 */
	const std::string summary () const;

private:
	double fitness = .0, meanSimilarity = .0, distancePenalty = .0;
	ConformerMatcher & conformerMatcher;
	DifgapeChromosome(const DifgapeChromosome & other) = delete;
	DifgapeChromosome & operator =(const DifgapeChromosome & other) = delete;
}
;

} /* namespace Difgape */

#endif /* DIFGAPECHROMOSOME_H_ */
