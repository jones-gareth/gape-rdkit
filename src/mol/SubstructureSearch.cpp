/*
 * SubstructureSearch.cpp
 *
 *  Created on: Oct 23, 2015
 *      Author: gjones
 */

#include "SubstructureSearch.h"

namespace GarethMol {

using namespace std;
using namespace GarethUtil;

SubstructureSearch::SubstructureSearch(const QueryMolecule & query, const Molecule & target) :
		queryMolecule(query), targetMolecule(target) {
	callbackFunction = [this](const vector<size_t> & queryIdsToTargetIds) {
		defaultCallbackFunction(queryIdsToTargetIds);
	};
	for (auto & queryAtom : queryMolecule.getAtoms()) {
		queryAtom->getSmartsAst()->initialize(targetMolecule);
	}
}

size_t SubstructureSearch::compare() {
	// Find all isomorphisms
	Ullmann<SubstructureSearch> ullmann(*this);
	auto nIsomorphisms = ullmann.doUllman();

	REPORT(Reporter::TRACE) << "found " << nIsomorphisms
			<< " structure isomorphisms";
	return nIsomorphisms;
}

size_t SubstructureSearch::getSizeOfTarget() const {
	return targetMolecule.nAtoms();
}

size_t SubstructureSearch::getSizeOfQuery() const {
	return queryMolecule.getAtoms().size();
}

bool SubstructureSearch::isAdjacentInTarget(const size_t targetNodeNo1,
		const size_t targetNodeNo2) const {
	const auto & bond = targetMolecule.getBond(targetNodeNo1, targetNodeNo2);
	return bond != nullptr;
}

bool SubstructureSearch::isAdjacentInQuery(const size_t queryNodeNo1,
		const size_t queryNodeNo2) const {
	const auto bond = queryMolecule.getBond(queryNodeNo1, queryNodeNo2);
	return bond != nullptr;
}

bool SubstructureSearch::sameNodeType(const size_t queryNodeNo,
		const size_t targetNodeNo) const {
	// check for fixed atoms
	auto iter = fixedAtomMapping.find(queryNodeNo);
	if (iter != fixedAtomMapping.end()) {
		auto testTargetNodeNo = iter->second;
		if (targetNodeNo != testTargetNodeNo) {
			return false;
		}
	}

	auto & queryAtom = queryMolecule.getAtoms().at(queryNodeNo);
	auto & targetAtom = targetMolecule.getAtom(targetNodeNo);
	auto match = queryAtom->matchAtom(targetAtom, targetMolecule);
	return match;
}

bool SubstructureSearch::sameEdgeType(const size_t queryNodeNo1,
		const size_t queryNodeNo2, const size_t targetNodeNo1,
		const size_t targetNodeNo2) const {
	const auto queryBond = queryMolecule.getBond(queryNodeNo1, queryNodeNo2);
	assert (queryBond != nullptr);
	const auto & targetBond = targetMolecule.getBond(targetNodeNo1,targetNodeNo2);
	assert (targetBond != nullptr);
	auto match = queryBond->matchBond(*targetBond, targetMolecule);
	return match;
}

void SubstructureSearch::callback(
		const vector<size_t> & queryIdsToTargetIds) const {
	callbackFunction(queryIdsToTargetIds);
}

void SubstructureSearch::defaultCallbackFunction(
		const vector<size_t> & queryIdsToTargetIds) const {
	REPORT(Reporter::TRACE) << "Found isomorphism " << endl;
	for (size_t queryId = 0; queryId < queryIdsToTargetIds.size(); queryId++) {
		auto targetId = queryIdsToTargetIds.at(queryId);
		auto & targetAtom = targetMolecule.getAtom(targetId);
		REPORT(Reporter::TRACE) << "Query atom " << queryId << " maps to target atom "
				<< targetAtom.info();
	}
}

} /* namespace GarethMol */
