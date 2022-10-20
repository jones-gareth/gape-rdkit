/*
 * ConformerPair.h
 *
 *  Created on: Jul 26, 2014
 *      Author: Gareth Jones
 */

#ifndef CONFORMERPAIR_H_
#define CONFORMERPAIR_H_

namespace Difgape {

class ConformerPair {
public:
	ConformerPair(int s1, int s2, int c1, int c2);
	virtual ~ConformerPair();
	inline bool operator< (const ConformerPair & rhs) const;

private:
	const int structure1, structure2, conformer1, conformer2;
	ConformerPair(const ConformerPair & other) = delete;
	ConformerPair & operator = (const ConformerPair & other) = delete;
};

class ConformerTriplet {

	ConformerTriplet(int s1, int s2, int s3, int c1, int c2, int c3);
	virtual ~ConformerTriplet();
	inline bool operator<(const ConformerTriplet & rhs) const;
private :
	const int structure1, structure2, structure3, conformer1, conformer2, conformer3;
	ConformerTriplet(const ConformerTriplet & ) = delete;
	ConformerTriplet operator = (const ConformerTriplet &) = delete;
};

class ScoreBounds {
	public:
	ScoreBounds (double min, double max);
	virtual ~ScoreBounds() {;}

	const double getMaxScore() const {
		return maxScore;
	}

	const double getMinScore() const {
		return minScore;
	}

private:
	const double minScore, maxScore;


};

} /* namespace Difgape */

#endif /* CONFORMERPAIR_H_ */
