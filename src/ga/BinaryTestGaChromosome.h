//
//  Copyright (C) 2020 Gareth Jones, Glysade LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

/**
 * A class to represent a Chromosome used in the binary f6 problem
 */

#ifndef BINARYTESTGACHROMOSOME_H_
#define BINARYTESTGACHROMOSOME_H_

#include "StringChromosome.h"

namespace Gape {

class BinaryTestGaChromosome: public BinaryStringChromosome {
private:
	double fitness = .0, xVal = .0, yVal = .0;
public:
    BinaryTestGaChromosome(const BinaryTestGaChromosome & other) = delete;
    BinaryTestGaChromosome & operator =(const BinaryTestGaChromosome & other) = delete;

	BinaryTestGaChromosome(RandomUtil & rng_,
			BinaryStringChromosomePolicy & chromosomePolicy_) :
			BinaryStringChromosome(44, rng_, chromosomePolicy_) {
	}

	~BinaryTestGaChromosome() = default;

	double getFitness() const {
		return fitness;
	}

	void clean() {
		fitness = .0;
		xVal = .0;
		yVal = .0;
	}

	bool isOk() {
		return true;
	}

	double score();

	double getXVal() const {
		return xVal;
	}

	double getYVal() const {
		return yVal;
	}

	std::string info() const;

	static bool sameNiche(const BinaryTestGaChromosome &other) {return false;}
};

}

#endif /* BINARYTESTGACHROMOSOME_H_ */
