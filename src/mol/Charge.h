/*
 * Charge.h
 *
 *  Created on: Jan 21, 2016
 *      Author: gjones
 */

#ifndef SRC_MOL_CHARGE_H_
#define SRC_MOL_CHARGE_H_

#include <string>
#include <vector>
#include <memory>
#include <SmartsParser.h>

namespace GarethMol {

using namespace std;

/**
 * A class to store charge patterns
 */
class ChargePattern {
	friend class Charge;

public:
	ChargePattern(SmartsParser & parser, string & patt, int fc, double pc,
			bool m) :
			pattern(patt), formalCharge(fc), partialCharge(pc), multiple(m) {
		query = parser.parseSmarts(patt);
	}

	virtual ~ChargePattern() {
	}

	ChargePattern(const ChargePattern & rhs) = delete;
	ChargePattern & operator =(const ChargePattern & rhs) = delete;
	ChargePattern(ChargePattern && rhs) = delete;
	ChargePattern & operator =(ChargePattern && rhs) = delete;

private:

	const string pattern;
	const int formalCharge;
	const double partialCharge;
	const bool multiple;
	unique_ptr<QueryMolecule> query;
};

/**
 * A class to handle the identification of charge patterns in a molecule.
 *
 */
class Charge {
public:
	Charge() {
	}

	virtual ~Charge() {
	}

	Charge(const Charge & rhs) = delete;
	Charge & operator =(const Charge & rhs) = delete;
	Charge(Charge && rhs) = delete;
	Charge & operator =(Charge && rhs) = delete;

	/**
	 * Apply charge patterns to a molecule
	 *
	 * @param molecule
	 */
	void chargeMolecule(Molecule & molecule);
private:

	/**
	 * Handle a matched charge pattern
	 *
	 * @param pattern
	 * @param molecule
	 * @param queryIdsToTargetIds
	 */
	void setChargeGroup(const ChargePattern & pattern, Molecule & molecule,
			const std::vector<size_t> & queryIdsToTargetIds) const;

	static std::vector<unique_ptr<ChargePattern>> & chargePatterns();
};

} /* namespace GarethMol */

#endif /* SRC_MOL_CHARGE_H_ */
