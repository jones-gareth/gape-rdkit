/*
 * RocsSdfMolecule.cpp
 *
 *  Created on: May 14, 2014
 *      Author: Gareth Jones
 */

#include <regex>
#include "RocsSdfMolecule.h"

#include "../util/Util.h"
#include "Reporter.h"
#include "RocsSdfParser.h"

namespace Difgape {

using namespace std;
using namespace GarethMol;
using namespace GarethUtil;

static regex conformerMatchPattern(R"(^(.*)_(\d+)\z)");

/**
 * Determines conformer number from molecule name
 *
 * @param molecule
 * @return
 */
static const int moleculeToConformerNo(const Molecule::MoleculePtr & molecule) {
	const string & name = molecule->getName();
	smatch matches;
	if (regex_search(name, matches, conformerMatchPattern)) {
		int conformerNo = convertString<int>(matches[2]);
		return conformerNo;
	} else {
		return 0;
	}
}

/**
 * Determines conformer name from molecule name
 * @param molecule
 * @return
 */
static const string moleculeToConformerName(
		const Molecule::MoleculePtr & molecule) {
	const string & name = molecule->getName();
	smatch matches;
	if (regex_search(name, matches, conformerMatchPattern)) {
		return matches[1];
	} else {
		return name;
	}
}

RocsSdfMolecule::RocsSdfMolecule(const Molecule::MoleculePtr & molecule) :
		structureName { moleculeToConformerName(molecule) }, conformerNo {
				moleculeToConformerNo(molecule) } {
	findFeatures(molecule);
	REPORT(Reporter::DETAIL) << "created ROCS sdf summary molecule "
			<< structureName << " conformer " << conformerNo;

}

void RocsSdfMolecule::findFeatures(const Molecule::MoleculePtr & molecule) {
	int nAtoms = molecule->getAtoms().size();
	featureCoordinates.resize(nAtoms, 4);
	featureNames.reserve(nAtoms);
	int no = 0;
	for (auto & atom : molecule->getAtoms()) {
		if (AtomType::isHeteroType(atom->getAtomType().getType())) {
			CoordVector c = molecule->getCoord(no);
			featureCoordinates.col(no) = c;
			featureNames.push_back(atom->info());
			no++;
		}
	}
	featureCoordinates.resize(no, 4);
}

bool RocsSdfMolecule::checkReference(RocsSdfMoleculePtr & other) const {
	const int noFeatures = this->noFeatures();
	if (noFeatures != other->noFeatures()) {
		return false;
	}
	for (int no = 0; no < noFeatures; no++) {
		if (getFeatureName(no) != other->getFeatureName(no)) {
			return false;
		}
	}
	return true;
}

void RocsSdfMolecule::setScore(Difgape::ScoreName scoreName, double score) {
	REPORT(Reporter::TRACE) << "Adding score field "
			<< RocsSdfParser::scoreNameToString(scoreName) << " value " << score
			<< " to molecule";
	assert(score >= 0);
	scores[scoreName] = score;
}

} /* namespace Difgape */
