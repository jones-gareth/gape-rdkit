/*
 * RocsSdfParser.h
 *
 *  Created on: May 14, 2014
 *      Author: Gareth Jones
 */

#ifndef ROCSSDFPARSER_H_
#define ROCSSDFPARSER_H_

#include <map>
#include <string>
#include <vector>

#include "RocsSdfMolecule.h"
#include "RocsSdfPair.h"

namespace Difgape {

using namespace std;

/**
 * Available score names
 */
enum class ScoreName {
	SHAPE_TANIMOTO,
	TVERSKY_Q,
	TVERSKY_D,
	SCALED_COLOR,
	COMBO_SCORE,
	OVERLAP,
	FIT_TVERSKY_COMBO,
	FIT_TVERSKY,
	FIT_COLOR_TVERSKY,
	REF_TVERSKY_COMBO,
	REF_TVERSKY,
	REF_COLOR_TVERSKY,
	SUB_TAN,
	First = SHAPE_TANIMOTO,
	Last = SUB_TAN
};

class RocsSdfParser {
public:
	explicit RocsSdfParser() {
	}
	;
	virtual ~RocsSdfParser() {
	}
	;

	/**
	 * Converts a score name to a string
	 *
	 * @param scoreName
	 * @return
	 */
	static const string scoreNameToString(const ScoreName & scoreName);

	/**
	 * converts a string to a score name
	 * @param label
	 * @return
	 */
	static const ScoreName stringToScoreName(const string & label);

	/**
	 * Returns the similarity score for a conformer pair
	 * @param scoreName
	 * @param compound1
	 * @param compound2
	 * @param conformer1
	 * @param conformer2
	 *
	 * @return
	 */
	const double getScore(ScoreName scoreName, int compound1, int compound2,
			int conformer1, int conformer2);

	/**
	 * Returns distance penalty score for a conformer triplet.
	 *
	 * @param compound1
	 * @param compound2
	 * @param compound3
	 * @param conformer1
	 * @param conformer2
	 * @param conformer3
	 * @return
	 */
	const double getPenaltyDistance(int compound1, int compound2, int compound3,
			int conformer1, int conformer2, int conformer3);

	/**
	 * Parses all files in a directory
	 *
	 * @param directory
	 */
	void readDirectory(const string & directory);

	const int getNoCompounds() const {
		return compoundNames.size();
	}

	const string & getCompoundName(int compoundNo) const {
		return compoundNames.at(compoundNo);
	}

	const int getConformerCount(int compoundNo) const {
		return conformerNumbers.at(compoundNo);
	}

	const string getConformerName(int compoundNo, int conformerNo) const {
		assert(conformerNo < conformerNumbers.at(compoundNo));
		return compoundNames.at(compoundNo) + "_" + to_string(conformerNo);
	}

	RocsSdfPairPtr getPair(string conformerName1, string conformerName2) const {
		return pairs.at(conformerName1).at(conformerName2);
	}

private:
	/**
	 * lookup of core structures. Index is structure name
	 */
	std::map<string, RocsSdfMolecule::RocsSdfMoleculePtr> referenceMolecules { };

	/**
	 * lookup up of pairs. Indexes are conformer names
	 */
	map<string, map<string, RocsSdfPairPtr>> pairs { };

	/**
	 * list of compound names, ordered by compound number
	 */
	vector<string> compoundNames { };

	/**
	 * stores total number of conformers by compound number
	 */
	vector<int> conformerNumbers { };

	/**
	 * Maps compound names to compound numbers
	 */
	map<string, int> compoundNumbers { };

	/**
	 * Reads and parses a ROCS results sdf file
	 * @param fileName
	 */
	void readFile(const string & fileName);

	/**
	 * Converts a regular molecule to a Rocs SDF molecule, performing checks and updating information in this class.
	 * @param mol
	 * @return
	 */
	RocsSdfMolecule::RocsSdfMoleculePtr toRocsSdfMolecule(
			const Molecule::MoleculePtr & mol);

};

} /* namespace Difgape */

#endif /* ROCSSDFPARSER_H_ */
