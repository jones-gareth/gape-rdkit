/*
 * SubstructureSearch.h
 *
 *  Created on: Oct 23, 2015
 *      Author: gjones
 */

#ifndef SRC_MOL_SUBSTRUCTURESEARCH_H_
#define SRC_MOL_SUBSTRUCTURESEARCH_H_

#include <memory>
#include "Molecule.h"
#include "QueryMolecule.h"
#include "../util/Ullmann.h"
#include "../util/Array2D.h"

namespace GarethMol {

using namespace std;
using namespace GarethUtil;

/**
 * Class to match a query molecule to a structure
 *
 * See Ullman.h
 */
class SubstructureSearch {
	friend Ullmann<SubstructureSearch> ;

public:
	SubstructureSearch(const QueryMolecule & query, const Molecule & target) ;

	virtual ~SubstructureSearch() {
	}

	SubstructureSearch(const SubstructureSearch & rhs) = delete;
	SubstructureSearch & operator =(const SubstructureSearch & rhs) = delete;
	SubstructureSearch(SubstructureSearch && rhs) = delete;
	SubstructureSearch & operator =(SubstructureSearch && rhs) = delete;

	/**
	 * Does the comparison and returns the number of subgraph isomorphisms found.
	 * The callback function will be called for each isomorphism found.
	 *
	 * @return
	 */
	size_t compare();

	/**
	 * Set the callback function required
	 *
	 * @param callbackFunction
	 */
	void setCallbackFunction(
			const function<void(const vector<size_t> & queryIdsToTargetIds)>& callbackFunction) {
		this->callbackFunction = callbackFunction;
	}

	void fixMapping(size_t queryAtomNo, size_t targetAtomNo) {
		fixedAtomMapping[queryAtomNo] =  targetAtomNo;
	}

	void clearFixedMappings() {
		fixedAtomMapping.clear();
	}
private:

	// These functions are required by the Ullman algorithm
	size_t getSizeOfTarget() const;
	size_t getSizeOfQuery() const;
	bool isAdjacentInTarget(const size_t targetNodeNo1,
			const size_t targetNodeNo2) const;
	bool isAdjacentInQuery(const size_t queryNodeNo1,
			const size_t queryNodeNo2) const;
	bool sameNodeType(const size_t queryNodeNo,
			const size_t targetNodeNo) const;
	bool sameEdgeType(const size_t queryNodeNo1, const size_t queryNodeNo2,
			const size_t targetNodeNo1, const size_t targetNodeNo2) const;
	void callback(const vector<size_t> & queryIdsToTargetIds) const;

	// this mapping fixes one or more query atoms to target atoms- their type still has to match.
	map<size_t, size_t> fixedAtomMapping {};

	const QueryMolecule & queryMolecule;
	const Molecule & targetMolecule;
	void defaultCallbackFunction(
			const vector<size_t> & queryIdsToTargetIds) const;
	function<void(const vector<size_t> & queryIdsToTargetIds)> callbackFunction;
};

} /* namespace GarethMol */

#endif /* SRC_MOL_SUBSTRUCTURESEARCH_H_ */
