/*
 * ConformerMatcher.h
 *
 *  Created on: Jun 11, 2014
 *      Author: Gareth Jones
 */

#ifndef CONFORMERMATCHER_H_
#define CONFORMERMATCHER_H_

#include "GaBase.h"
#include "GaOperation.h"
#include "RocsSdfParser.h"
#include "IntegerStringChromosomePolicy.h"
#include "LinkedPopLinearSel.h"
#include "DifgapeChromosome.h"

class DifgapeChromosome;

namespace Difgape {

using namespace GapeGa;

class ConformerMatcher;

using ConformerMatcherPopulation = LinkedPopLinearSel<DifgapeChromosome, ConformerMatcher>;

/**
 * Main class for running DIFGAPE
 */
class ConformerMatcher: public GaBase {
public:
	explicit ConformerMatcher(RocsSdfParser & parser);
	virtual ~ConformerMatcher() {
	}

	IntegerStringChromosomePolicy & getChromosomePolicy() {
		return chromosomePolicy;
	}

	RocsSdfParser & getRocsSdfParser() const {
		return rocsSdfParser;
	}

	ScoreName getScoreField() const {
		return scoreField;
	}

	const double getPenaltyWeight() const {
		return penaltyWeight;
	}

	/**
	 * Runs the GA
	 */
	void run();

	/**
	 * Create a difgape chromosome
	 *
	 * @return
	 */
	shared_ptr<DifgapeChromosome> createChromosome();

	const vector<shared_ptr<GaOperation<DifgapeChromosome> > > getOperations() const {
		return operations;
	}

private:
	RocsSdfParser & rocsSdfParser;
	IntegerStringChromosomePolicy chromosomePolicy;
	ScoreName scoreField = ScoreName::COMBO_SCORE;
	double penaltyWeight = 1.0;
	int nIslands = 5;
	double migrateWeight = 10;
	double crossoverWeight = 95;
	double mutateWeight = 95;
	int noIterations = 60000;
	int noRuns = 10;

	vector<shared_ptr<GaOperation<DifgapeChromosome> > > operations;
	unique_ptr<ConformerMatcherPopulation> population;

	/**
	 * Applies a mutation operator
	 *
	 * @param parents
	 * @param children
	 */
	static void mutate(const vector<shared_ptr<DifgapeChromosome> > & parents,
			vector<shared_ptr<DifgapeChromosome> > & children);

	/*
	 * Applies a crossover operator
	 */
	static void crossover(
			const vector<shared_ptr<DifgapeChromosome> > & parents,
			vector<shared_ptr<DifgapeChromosome> > & children);

	/**
	 * Creates GA operations
	 */
	void createOperations();

};

} /* namespace Difgape */

#endif /* CONFORMERMATCHER_H_ */
