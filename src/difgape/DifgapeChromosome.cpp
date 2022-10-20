/*
 * DifgapeChromosome.cpp
 *
 *  Created on: Jun 11, 2014
 *      Author: Gareth Jones
 */

#include "DifgapeChromosome.h"
#include "ConformerMatcher.h"
#include "boost/format.hpp"

namespace Difgape {

using namespace std;
using namespace GapeGa;

DifgapeChromosome::DifgapeChromosome(ConformerMatcher & cm_) :
		IntegerStringChromosome(cm_.getRocsSdfParser().getNoCompounds(),
				cm_.getRng(), cm_.getChromosomePolicy()), conformerMatcher(cm_) {
	;
}

double DifgapeChromosome::score() {

	double sum = 0;
	double nValues = 0;
	double nPairs = 0;
	RocsSdfParser & parser = conformerMatcher.getRocsSdfParser();
	ScoreName scoreField = conformerMatcher.getScoreField();

	for (int i = 0; i < length; i++) {
		const int conformer1 = getValue(i);
		if (conformer1 == -1) {
			continue;
		}
		nValues += 1.0;
		double score = .0;
		for (int j = i + 1; j < length; j++) {
			const int conformer2 = getValue(j);
			if (conformer2 == -1)
				continue;
			score += parser.getScore(scoreField, i, j, conformer1, conformer2);
			nPairs += 1.0;
		}
		sum += score;
	}

	meanSimilarity = nPairs > 0 ? sum / nPairs : 0;

	// check all triplets for distance penalties
	double penalty = 0, cnt = 0;
	for (int i = 0; i < length; i++) {
		const int conformer1 = getValue(i);
		if (conformer1 == -1)
			continue;
		for (int j = i + 1; j < length; j++) {
			const int conformer2 = getValue(j);
			if (conformer2 == -1)
				continue;
			for (int k = j + 1; k < length; k++) {
				const int conformer3 = getValue(k);
				if (conformer3 == -1)
					continue;
				cnt++;
				penalty += parser.getPenaltyDistance(i, j, k, conformer1,
						conformer2, conformer3);
			}
		}
	}
	distancePenalty = cnt > 0 ? penalty / cnt : 0;

	fitness = meanSimilarity
			- conformerMatcher.getPenaltyWeight() * distancePenalty;
	fitness *= nValues;
	return fitness;
}

const string DifgapeChromosome::info() const {
	boost::format format = boost::format(
			"Fit %6.3f sim %6.3f distance %6.3f : ") % fitness % meanSimilarity
			% distancePenalty;
	return format.str() + geneInfo();
}

const string DifgapeChromosome::summary() const {
	stringstream ss;
	RocsSdfParser & parser = conformerMatcher.getRocsSdfParser();
	ScoreName scoreField = conformerMatcher.getScoreField();

	ss<<"Fitness : " << getFitness() << endl;
	ss <<"Mean similarity : "<< meanSimilarity<< endl;
	ss<<"Distance penalty : "<<distancePenalty<<endl;
	ss<<"Mapped ";
	for (int i = 0; i < length; i++) {
		if (getValue(i)!=-1) {
			ss<<parser.getConformerName(i, getValue(i))<<" ";
		}
	}
	ss<<endl;
	ss<<"UnMapped ";
	for (int i = 0; i < length; i++) {
		if (getValue(i)==-1) {
			ss<<parser.getCompoundName(i)<<" ";
		}
	}
	ss<<endl;

	for (int i = 0; i < length; i++) {
		if (getValue(i) == -1)
			continue;
		for (int j = i + 1; j < length; j++) {
			if (getValue(j) == -1)
				continue;
			double score = parser.getScore(scoreField, i, j, getValue(i),
						getValue(j));

			ss << parser.getConformerName(i, getValue(i))<<"\t";
			ss << parser.getConformerName(j, getValue(j))<<"\t";
			ss<<score<<endl;

		}
	}

	return ss.str();
}

} /* namespace Difgape */
