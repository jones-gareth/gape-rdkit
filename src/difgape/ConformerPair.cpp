/*
 * ConformerPair.cpp
 *
 *  Created on: Jul 26, 2014
 *      Author: Gareth Jones
 */

#include "ConformerPair.h"
#include <cassert>

namespace Difgape {

ConformerPair::ConformerPair(int s1, int s2, int c1, int c2) :
		structure1(s1), structure2(s2), conformer1(c1), conformer2(c2) {
	assert(s1 < s2);
}

bool ConformerPair::operator< (
		const ConformerPair & rhs) const {
	if (this->structure1 != rhs.structure1)
		return this->structure1 < rhs.structure1;
	if (this->structure2 != rhs.structure2)
		return this->structure2 < rhs.structure2;
	if (this->conformer1 != rhs.conformer1)
		return this->conformer1 < rhs.conformer1;
	if (this->conformer2 != rhs.conformer2)
		return this->conformer2 < rhs.conformer2;
	return false;
}

ConformerPair::~ConformerPair() {
}

ConformerTriplet::ConformerTriplet(int s1, int s2, int s3, int c1, int c2,
		int c3) :
		structure1(s1), structure2(s2), structure3(s3), conformer1(c1), conformer2(
				c2), conformer3(c3) {
	assert(s1 < s2);
	assert(s2 < s3);
}

ConformerTriplet::~ConformerTriplet() {

}

bool ConformerTriplet::operator< (
		const ConformerTriplet & rhs) const {
	if (this->structure1 != rhs.structure1)
		return this->structure1 < rhs.structure1;
	if (this->structure2 != rhs.structure2)
		return this->structure2 < rhs.structure2;
	if (this->structure3 != rhs.structure3)
		return this->structure3 < rhs.structure3;
	if (this->conformer1 != rhs.conformer1)
		return this->conformer1 < rhs.conformer1;
	if (this->conformer2 != rhs.conformer2)
		return this->conformer2 < rhs.conformer2;
	if (this->conformer3 != rhs.conformer3)
		return this->conformer3 < rhs.conformer3;
	return false;


}

ScoreBounds::ScoreBounds(double min, double max) : minScore(min), maxScore(max) {
	assert(minScore<=maxScore);
}

}
/* namespace Difgape */
