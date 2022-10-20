/*
 * RocsSdfParser.cpp
 *
 *  Created on: May 14, 2014
 *      Author: Gareth Jones
 */

#include <boost/regex.hpp>
#include <boost/filesystem.hpp>

#include "RocsSdfParser.h"
#include "Util.h"
#include "Reporter.h"
#include "EnumIter.h"

namespace Difgape {

using namespace std;
using namespace GarethUtil;

static const string scoreNameLabels[] = { "ShapeTanimoto", "Tversky(q)",
		"Tversky(d)", "ScaledColor", "ComboScore", "Overlap", "FitTverskyCombo",
		"FitTversky", "FitColorTversky", "RefTverskyCombo", "RefTversky",
		"RefColorTversky", "SubTan" };

const string RocsSdfParser::scoreNameToString(const ScoreName & scoreName) {
	return scoreNameLabels[static_cast<int>(scoreName)];
}

const ScoreName RocsSdfParser::stringToScoreName(const string & label) {
	int nScoreNames = static_cast<int>(ScoreName::Last);
	for (int i = 0; i < nScoreNames; i++) {
		if (equalsIgnoreCase(scoreNameLabels[i], label)) {
			return static_cast<ScoreName>(i);
		}
	}
	throw invalid_argument("no score field with name " + label);
}

void RocsSdfParser::readDirectory(const string & directory) {
	REPORT(Reporter::INFO) << "Reading directory " << directory;
	if (!boost::filesystem::exists(directory)) {
		throw invalid_argument(directory + " does not exist!");
	}
	boost::filesystem::directory_iterator end_itr; // default construction yields past-the-end
	for (boost::filesystem::directory_iterator itr(directory); itr != end_itr;
			++itr) {
		if (boost::filesystem::is_regular_file(itr->status())) {
			string fileName { itr->path().string() };
			string test = toLowerCase(fileName);
			if ((endsWith(test, "sdf") || endsWith(test, "sdf.gz"))
					&& test.find("hits") != string::npos) {
				readFile(fileName);
			}
		}
	}

	for (int i = 0; i < getNoCompounds(); i++) {
		REPORT(Reporter::INFO) << "Structure " << (i + 1) << ": "
				<< compoundNames.at(i) << " # conformers "
				<< conformerNumbers.at(i);
	}
	REPORT(Reporter::INFO) << "Finished reading directory";
}

void RocsSdfParser::readFile(const string & fileName) {
	REPORT(Reporter::INFO) << "Reading file " << fileName;
	vector<Molecule::MoleculePtr> mols = Molecule::readMoleculesFromFile(
			fileName);

	Molecule::MoleculePtr & firstMolecule = mols.at(0);
	// note the first molecule name may not contain the conformer name
	string queryStructureName = firstMolecule->getName();

	// first we need to find the query molecule. This is normally the first
	// match, but if we have identical shapes in the input structures it
	// may be further down- check the first 10 for matching structure names.
	RocsSdfMolecule::RocsSdfMoleculePtr queryMolecule = nullptr;
	for (size_t i = 1; i < 11; i++) {
		if (mols.size() == i)
			break;
		RocsSdfMolecule::RocsSdfMoleculePtr targetMolecule { toRocsSdfMolecule(
				mols.at(i)) };
		string targetMoleculeName = targetMolecule->getStructureName();
		if (targetMoleculeName == queryStructureName) {
			REPORT(Reporter::DEBUG) << "Entry " << i << " matches query";
			queryMolecule = targetMolecule;
			break;
		}
	}

	if (queryMolecule == nullptr)
		throw runtime_error("can't read query molecule");
	string queryName = queryMolecule->getName();
	REPORT(Reporter::INFO) << "Query is " << queryName;
	RocsSdfMolecule::RocsSdfMoleculePtr queryReference = referenceMolecules.at(
			queryMolecule->getStructureName());
	if (!isKeyPresent(queryName, pairs)) {
		pairs[queryName];
	}

	// process file, skipping first entry
	for (size_t i = 1; i < mols.size(); i++) {
		RocsSdfMolecule::RocsSdfMoleculePtr targetMolecule { toRocsSdfMolecule(
				mols.at(i)) };
		string targetName = targetMolecule->getName();
		// check target isn't query
		if (targetName != queryName) {

			RocsSdfMolecule::RocsSdfMoleculePtr targetReference =
					referenceMolecules.at(targetMolecule->getStructureName());

			// create pair structure
			RocsSdfPairPtr pair { new RocsSdfPair(queryMolecule, targetMolecule,
					fileName) };
			REPORT(Reporter::DEBUG) << "reading pair " << queryName << ","
					<< targetName;

			// add pair to lookup table
			if (isKeyPresent(targetName, pairs.at(queryName))) {
				throw runtime_error(
						"pair for " + queryName + " and " + targetName
								+ " already defined");
			}
			pairs.at(queryName)[targetName] = pair;
		}
	}
	RocsSdfMolecule::RocsSdfMoleculePtr rocsSdfMol;

}

RocsSdfMolecule::RocsSdfMoleculePtr RocsSdfParser::toRocsSdfMolecule(
		const Molecule::MoleculePtr & mol) {
	if (mol == nullptr)
		return nullptr;
	RocsSdfMolecule::RocsSdfMoleculePtr molecule { new RocsSdfMolecule(mol) };

	// add scores - this should maybe be done in RocsSdfMolecule constructor
	for (ScoreName scoreName : EnumIter<ScoreName>()) {
		string strVal = RocsSdfParser::scoreNameToString(scoreName);
		REPORT(Reporter::TRACE) << "Checking score field " << strVal;
		for (string str : { strVal, "ROCS_" + strVal }) {
			if (mol->hasSdfField(str)) {
				molecule->setScore(scoreName,
						stod(mol->getSdfValues(str).at(0)));
			}
		}
	}

	// add to reference molecule list or make sure it's consistent with current
	// reference molecule
	string name = molecule->getStructureName();
	auto pair = referenceMolecules.find(name);
	if (pair == referenceMolecules.end()) {
		referenceMolecules[name] = molecule;
	} else {
		RocsSdfMolecule::RocsSdfMoleculePtr refMolecule = pair->second;
		if (!molecule->checkReference(refMolecule)) {
			throw runtime_error("Reference error for " + molecule->getName());
		}
	}

	// this is the number of conformers available (rather that the conformer
	// number)!
	int conformerNo = molecule->getConformerNo() + 1;
	if (!contains(compoundNames, name)) {
		// haven't seen this structure before. Update compound names list
		// and conformer numbers list
		int moleculeNo = compoundNames.size();
		compoundNames.push_back(name);
		compoundNumbers[name] = moleculeNo;
		conformerNumbers.push_back(conformerNo);
	} else {
		// we know about this compound. Update conformer numbers if this is
		// the largest conformer number we've seen.
		int moleculeNo = compoundNumbers.at(name);
		assert(compoundNames.at(moleculeNo) == name);
		if (conformerNumbers.at(moleculeNo) < conformerNo)
			conformerNumbers.at(moleculeNo) = conformerNo;
	}

	return molecule;
}

const double RocsSdfParser::getScore(ScoreName scoreName, int compound1,
		int compound2, int conformer1, int conformer2) {
	string queryName = getConformerName(compound1, conformer1);
	string targetName = getConformerName(compound2, conformer2);
	return pairs.at(queryName).at(targetName)->getScore(scoreName);
}

const double RocsSdfParser::getPenaltyDistance(int compound1, int compound2,
		int compound3, int conformer1, int conformer2, int conformer3) {
	// get full structure names
	string structure1 = getConformerName(compound1, conformer1);
	string structure2 = getConformerName(compound2, conformer2);
	string structure3 = getConformerName(compound3, conformer3);

	int nFeatures1 =
			referenceMolecules.at(compoundNames.at(compound1))->getNFeatures();
	int nFeatures2 =
			referenceMolecules.at(compoundNames.at(compound2))->getNFeatures();
	int nFeatures3 =
			referenceMolecules.at(compoundNames.at(compound3))->getNFeatures();

	RocsSdfPairPtr pair12 = getPair(structure1, structure2);
	RocsSdfPairPtr pair23 = getPair(structure2, structure3);
	RocsSdfPairPtr pair31 = getPair(structure3, structure1);

	double penalty = 0;
	double cnt = 0;

	for (int i = 0; i < nFeatures1; i++) {
		for (int j = 0; j < nFeatures2; j++) {
			for (int k = 0; k < nFeatures3; k++) {

				double distance12 = pair12->getDistance(i, j);
				double distance23 = pair23->getDistance(j, k);
				double distance31 = pair31->getDistance(k, i);

				double max = 0, sum = 0;
				if (distance12 >= distance23 && distance12 >= distance31) {
					max = distance12;
					sum = distance23 + distance31;
				} else if (distance23 >= distance12
						&& distance23 >= distance31) {
					max = distance23;
					sum = distance12 + distance31;
				} else if (distance31 >= distance12
						&& distance31 >= distance23) {
					max = distance31;
					sum = distance12 + distance23;
				} else
					assert(false);

				assert(
						abs(distance12 + distance23 + distance31 - max - sum)
								< 1e-10);
				cnt++;
				if (max > sum) {
					double diff = max - sum;
					penalty += diff * diff;
				}
			}
		}
	}

	return penalty / cnt;
}

} /* namespace Difgape */
