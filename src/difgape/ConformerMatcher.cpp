/*
 * ConformerMatcher.cpp
 *
 *  Created on: Jun 11, 2014
 *      Author: Gareth Jones
 */

#include "ConformerMatcher.h"
#include "../util/Reporter.h"

namespace Difgape {

using namespace std;
using namespace GapeGa;

ConformerMatcher::ConformerMatcher(RocsSdfParser & parser) :
		rocsSdfParser(parser), chromosomePolicy(getRng(),
				parser.getNoCompounds()) {
	setPopsize(500);
	setSelectionPressure(1.0001);

	chromosomePolicy.setAllowNulls(true);
	for (int i = 0; i < parser.getNoCompounds(); i++) {
		chromosomePolicy.setMax(i, parser.getConformerCount(i));
	}

}

void ConformerMatcher::mutate(
		const vector<shared_ptr<DifgapeChromosome> > & parents,
		vector<shared_ptr<DifgapeChromosome> > & children) {
	shared_ptr<DifgapeChromosome> parent = parents[0];
	shared_ptr<DifgapeChromosome> child = children[0];

	child->copyGene(*parent);
	child->mutate();
}

void ConformerMatcher::crossover(
		const vector<shared_ptr<DifgapeChromosome> > & parents,
		vector<shared_ptr<DifgapeChromosome> > & children) {

	shared_ptr<DifgapeChromosome> parent1 = parents[0];
	shared_ptr<DifgapeChromosome> child1 = children[0];
	shared_ptr<DifgapeChromosome> parent2 = parents[1];
	shared_ptr<DifgapeChromosome> child2 = children[1];

	parent1->onePointCrossover(*parent2, *child1, *child2);
}

void ConformerMatcher::createOperations() {
	shared_ptr<GaOperation<DifgapeChromosome> > mutationOperation(
			new GaOperation<DifgapeChromosome>(1, 1, mutateWeight,
					&ConformerMatcher::mutate));
	shared_ptr<GaOperation<DifgapeChromosome> > crossoverOperation(
			new GaOperation<DifgapeChromosome>(2, 2, crossoverWeight,
					&ConformerMatcher::crossover));
	operations.reserve(2);
	operations.push_back(mutationOperation);
	operations.push_back(crossoverOperation);
}

shared_ptr<DifgapeChromosome> ConformerMatcher::createChromosome() {
	shared_ptr<DifgapeChromosome> chromosome(new DifgapeChromosome(*this));
	REPORT(Reporter::TRACE) << "Created chromosome at address "
			<< chromosome.get();
	return chromosome;
}

void ConformerMatcher::run() {
	createOperations();
	population = unique_ptr<ConformerMatcherPopulation>(
			new ConformerMatcherPopulation(*this));
	population->create();
	REPORT(Reporter::INFO) << population->info();
	REPORT(Reporter::DETAIL) << population->populationInfo();

	int nOps = 0;
	while (nOps < noIterations) {
		population->iterate();
		nOps++;
		if (nOps % 1000 == 0) {
			REPORT(Reporter::INFO) << population->info();
		}
	}
	const shared_ptr<DifgapeChromosome> best = population->getBest();
	REPORT(Reporter::INFO) << "Best solution " << best->info() << endl
			<< best->summary();
	REPORT(Reporter::DETAIL) << population->populationInfo();
}

} /* namespace Difgape */
