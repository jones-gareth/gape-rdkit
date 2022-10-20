/*
 * MolecularRms.cpp
 *
 *  Created on: Sep 10, 2015
 *      Author: gjones
 */

#include <boost/optional.hpp>
#include "Eigen/Dense"
#include "MolecularRms.h"
#include "../util/LeastSquaresFit.h"

namespace GarethMol {

using namespace std;
using namespace GarethUtil;

MolecularRms::~MolecularRms() {
}

void MolecularRms::determineIsomorphisms() {
	molComparer = make_unique<MolComparer>(queryMolecule, targetMolecule);
	molComparer->setHeavyAtomOnly(true);
	molComparer->setSubgraph(false);

	auto callbackFunction =
			[this] (const std::vector<size_t> & queryIdsToTargetIds) {
				addIsomorphism (queryIdsToTargetIds);
			};
	isomorphisms.clear();
	molComparer->setCallbackFunction(callbackFunction);
	molComparer->compare();
	if (isomorphisms.size() == 0) {
		REPORT(Reporter::FATAL) << "No isomorphisms found!";
		throw runtime_error("No isomorphisms found!");
	}
}

double MolecularRms::determineRms(const CoordMatrix & queryCoords,
		const CoordMatrix & targetCoords) {

	this->queryCoords = queryCoords;
	this->targetCoords = targetCoords;

	minRmsd = numeric_limits<double>::max();

	evaluateIsomorphisms();

	REPORT(Reporter::DETAIL) << "Rmsd values "
			<< collectionToString(isomorphismRmsds);
	return minRmsd;
}

void MolecularRms::evaluateIsomorphisms() {
	minRmsd = numeric_limits<double>::max();
	isomorphismRmsds.clear();

	const auto & queryAtoms = molComparer->getQueryAtoms();
	const auto & targetAtoms = molComparer->getTargetAtoms();

	for (const auto & queryIdsToTargetIds : isomorphisms) {
		auto rmsd = leastSquaresAtomsFit(queryIdsToTargetIds, queryAtoms,
				targetAtoms, queryCoords, targetCoords, doLeastSquaresFit);
		if (rmsd < minRmsd) {
			minRmsd = rmsd;
		}
		isomorphismRmsds.push_back(rmsd);
	}
}

unique_ptr<MolecularRms> MolecularRms::createSelfComparer(
		const Molecule & molecule) {
	auto molecularRms = make_unique<MolecularRms>(molecule, molecule);
	return molecularRms;
}

double MolecularRms::leastSquaresAtomsFit(
		const std::vector<size_t> & queryIdsToTargetIds,
		const std::vector<const Atom *> & queryAtoms,
		const std::vector<const Atom *> & targetAtoms,
		const CoordMatrix & queryCoords, const CoordMatrix & targetCoords,
		bool doLeastSquaresFit) {
	size_t size = queryIdsToTargetIds.size();
	CoordMatrix mappedQueryCoords(4, size), mappedTargetCoords(4, size);

	for (size_t queryId = 0; queryId < size; queryId++) {
		auto targetId = queryIdsToTargetIds.at(queryId);
		auto queryAtomNo = queryAtoms.at(queryId)->getAtomNo();
		auto targetAtomNo = targetAtoms.at(targetId)->getAtomNo();
		REPORT(Reporter::TRACE) << "query " << queryAtomNo
				<< " maps to target " << targetAtomNo;
		mappedQueryCoords.col(queryId) = queryCoords.col(queryAtomNo);
		mappedTargetCoords.col(queryId) = targetCoords.col(targetAtomNo);
	}

	REPORT(Reporter::TRACE) << "query coords " << endl << queryCoords;
	REPORT(Reporter::TRACE) << "mapped query coords " << endl
			<< mappedQueryCoords;
	REPORT(Reporter::TRACE) << "target coords " << endl << targetCoords;
	REPORT(Reporter::TRACE) << "mapped target coords " << endl
			<< mappedTargetCoords;

	double rmsd;
	if (doLeastSquaresFit) {
		auto matrix = leastSquaresFit(mappedQueryCoords, mappedTargetCoords,
				boost::none);

		auto & transformedQueryCoords = matrix * mappedQueryCoords;
		REPORT(Reporter::TRACE) << "transformed query coords " << endl
				<< transformedQueryCoords;
		rmsd = rmsDistance(transformedQueryCoords, mappedTargetCoords);
	} else {
		rmsd = rmsDistance(mappedQueryCoords, mappedTargetCoords);
	}
	return rmsd;
}

void MolecularRms::addIsomorphism(const std::vector<size_t> & queryIdsToTargetIds) {
	isomorphisms.push_back(queryIdsToTargetIds);
}


} /* namespace GarethUtil */
