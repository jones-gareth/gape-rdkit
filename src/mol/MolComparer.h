/*
 * MolComparer.h
 *
 *  Created on: Sep 3, 2015
 *      Author: gjones
 */

#ifndef SRC_MOL_MOLCOMPARER_H_
#define SRC_MOL_MOLCOMPARER_H_

#include <memory>
#include "Molecule.h"
#include "../util/Ullmann.h"

using namespace std;
using namespace GarethUtil;

namespace GarethMol {

/**
 * A class which uses the Ullman algorithm to evaluate all isomorphisms between
 * a query molecule and a target molecule.
 *
 * Can match on heavy atoms only, or can do subgraph isomorphism
 *
 * This class has a default callback function which just prints out the mapping found.
 *
 * The class implements the policy required by the Ullman template class.
 */
class MolComparer {
	friend Ullmann<MolComparer> ;

public:
	void defaultCallbackFunction(const std::vector<size_t> & queryIdsToTargetIds);

	MolComparer(const Molecule & queryMol, const Molecule & targetMol) :
			queryMolecule(queryMol), targetMolecule(targetMol) {
		callbackFunction = [this](const std::vector<size_t> & queryIdsToTargetIds) {
			defaultCallbackFunction(queryIdsToTargetIds);
		};
	}

	virtual ~MolComparer() {
	}

	MolComparer(const MolComparer & rhs) = delete;
	MolComparer & operator =(const MolComparer & rhs) = delete;
	MolComparer(MolComparer && rhs) = delete;
	MolComparer & operator =(MolComparer && rhs) = delete;

	/**
	 * Does the comparison and returns the number of isomorphisms found.
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
			const function<void(const std::vector<size_t> & queryIdsToTargetIds)>& callbackFunction) {
		this->callbackFunction = callbackFunction;
	}

	bool isHeavyAtomOnly() const {
		return heavyAtomOnly;
	}

	void setHeavyAtomOnly(bool heavyAtomOnly) {
		this->heavyAtomOnly = heavyAtomOnly;
	}

	bool isSubgraph() const {
		return subgraph;
	}

	void setSubgraph(bool subgraph) {
		this->subgraph = subgraph;
	}

	/**
	 * Create and return a comparer that is suitable for comparing a molecule with itself.
	 *
	 * @param molecule
	 * @return
	 */
	static unique_ptr<MolComparer> createIsomorphismComparer(
			const Molecule &molecule);

	const std::vector<const Atom *>& getQueryAtoms() const {
		return queryAtoms;
	}

	const std::vector<const Atom *>& getTargetAtoms() const {
		return targetAtoms;
	}

	size_t getIsomorphisms() const {
		return nIsomorphisms;
	}

	bool isMatchElementalTypes() const {
		return matchElementalTypes;
	}

	void setMatchElementalTypes(bool matchElementalTypes) {
		this->matchElementalTypes = matchElementalTypes;
	}

    bool isMatchTypes() const {
        return matchTypes;
    }

    void setMatchTypes(bool matchTypes = true) {
        this->matchTypes = matchTypes;
    }

private:
	bool heavyAtomOnly = false;
	bool subgraph = true;
	bool matchElementalTypes = true;
	bool matchTypes = true;

	// These functions are required by the Ullman algorithm
	size_t getSizeOfTarget();
	size_t getSizeOfQuery();
	bool isAdjacentInTarget(const size_t targetNodeNo1,
			const size_t targetNodeNo2);
	bool isAdjacentInQuery(const size_t queryNodeNo1,
			const size_t queryNodeNo2);
	bool sameNodeType(const size_t queryNodeNo,
			const size_t targetNodeNo) const;
	bool sameEdgeType(const size_t queryNodeNo1,
			const size_t queryNodeNo2, const size_t targetNodeNo1,
			const size_t targetNodeNo2) const;
	void callback(const std::vector<size_t> & queryIdsToTargetIds);

	std::vector<const Atom *> findAtoms(const Molecule & molecule) const;

	const Molecule & queryMolecule, &targetMolecule;
	size_t nIsomorphisms = 0;
	std::vector<const Atom *> queryAtoms, targetAtoms;
	std::vector<size_t> queryAtomConnections, targetAtomConnections;
	unique_ptr<Array2D<bool>> queryAdjacencyMatrix = nullptr,
			targetAdjacencyMatrix = nullptr;
	function<void(const std::vector<size_t> & queryIdsToTargetIds)> callbackFunction;

};

} /* namespace GarethMol */

#endif /* SRC_MOL_MOLCOMPARER_H_ */
