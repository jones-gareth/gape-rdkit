/*
 * Ullman.h
 *
 *  Created on: Sep 1, 2015
 *      Author: gjones
 */

#ifndef SRC_UTIL_ULLMAN_H_
#define SRC_UTIL_ULLMAN_H_

#include <cassert>
#include <vector>
#include <map>
#include <set>
#include "Reporter.h"

namespace GarethUtil {

/**
 * A template class to perform the Ullman subgraph isomorphism
 *
 * In order to use you need to instantiate using a policy class (member info) with the following methods:
 *
 * size_t getSizeOfTarget() const;
 * size_t getSizeOfQuery() const;
 * bool isAdjacentInTarget(const size_t targetNodeNo1, const size_t targetNodeNo2) const;
 * bool isAdjacentInQuery(const size_t queryNodeNo1, const size_t queryNodeNo2) const;
 * callback(const & std::vector<size_t> queryIdsToTargetIds)
 * bool sameNodeType(const size_t queryNodeNo, const size_t targetNodeNo) const;
 * bool sameEdgeType(const size_t queryNodeNo1, const size_t queryNodeNo2, const size_t targetNodeNo1,
 *                   const size_t targetNodeNo2) const;
 *
 * See http://stackoverflow.com/questions/13537716/how-to-partially-compare-two-graphs
 * for pretty decent Python implementation of the Ullman algorithm
 *
 */
template<typename UllmanInfo>
class Ullmann {
public:
	Ullmann(UllmanInfo & info_) :
			info(info_), querySize(info.getSizeOfQuery()), targetSize(
					info.getSizeOfTarget()), queryAdjacencies(querySize), targetAdjacencies(
					targetSize) {
	}

	virtual ~Ullmann() {
	}

	Ullmann(const Ullmann & rhs) = delete;
	Ullmann & operator =(const Ullmann & rhs) = delete;
	Ullmann(Ullmann && rhs) = delete;
	Ullmann & operator =(Ullmann && rhs) = delete;

	size_t doUllman() {
		nIsomorphisms = 0;
		search();
		return nIsomorphisms;
	}
private:
	using ListOfSets = std::vector<set<size_t>>;
	UllmanInfo & info;
	int nIsomorphisms = 0;
	const size_t querySize;
	const size_t targetSize;
	ListOfSets queryAdjacencies;
	ListOfSets targetAdjacencies;
	bool findAllIsoMorphisms = true;

	/**
	 * Initialize data structures and call depth first search
	 */
	void search() {
		createAdjacencies(true);
		createAdjacencies(false);
		std::vector<size_t> mapping;
		mapping.reserve(querySize);
		ListOfSets allowedMappings(querySize);
		// create initial allowed correspondances
		for (size_t queryNode = 0; queryNode < querySize; queryNode++) {
			for (size_t targetNode = 0; targetNode < targetSize; targetNode++) {
				if (info.sameNodeType(queryNode, targetNode)) {
					allowedMappings.at(queryNode).insert(targetNode);
				}
			}
		}

		depthFirst(mapping, allowedMappings);
	}

	/**
	 * Does the depth first search
	 *
	 * @param mapping
	 * @param allowedMappings
	 * @return
	 */
	bool depthFirst(std::vector<size_t> & mapping, ListOfSets & allowedMappings) {
		updateAllowedMappings(allowedMappings);

		auto size = mapping.size();
		// check that query matches target so far
		for (size_t i=0; i < size; i++) {
			for (size_t j : queryAdjacencies.at(i)) {
				if (j < size) {
					// nodes i and j are adjacent in the query graph
					size_t x = mapping.at(i);
					size_t y = mapping.at(j);
					// backtrack if target node pairs are not adjacent
					if (!targetAdjacencies.at(x).count(y))
						return false;
					// For labeled graphs check that query and target edges match
					if (!info.sameEdgeType(i, j, x, y))
						return false;
				}
			}
		}

		if (size == querySize) {
			// found a subgraph isomorphism
			nIsomorphisms++;
			info.callback(mapping);
			return !findAllIsoMorphisms;
		}

		// expand the search to the next depth/query node

		auto & allowedTargetNodes = allowedMappings.at(size);
		REPORT(Reporter::TRACE) << "allowed Target nodes: "
				<< collectionToString(allowedTargetNodes);
		auto allowedMappingsCopy = allowedMappings.at(size);
		for (auto nextTargetNode : allowedMappingsCopy) {
			REPORT(Reporter::TRACE) << "Checking node " << nextTargetNode;
			if (!contains(mapping, nextTargetNode)) {
				mapping.push_back(nextTargetNode);
				// Create a new set of possible assignments, where graph node j is the only
				// possibility for the assignment of subgraph node i.
				// newAllowedMappings *must* be a copy of allowedMappings
				auto newAllowedMappings = allowedMappings;
				newAllowedMappings.at(size) = {nextTargetNode};
				if (depthFirst(mapping, newAllowedMappings)) {
					return true;
				}
				mapping.pop_back();
			}
			REPORT(Reporter::TRACE) << "allowed Target nodes before erase: "
					<< collectionToString(allowedTargetNodes);
			allowedTargetNodes.erase(nextTargetNode);
			REPORT(Reporter::TRACE) << "allowed Target nodes after erase: "
					<< collectionToString(allowedTargetNodes);
			updateAllowedMappings(allowedMappings);
		}

		return false;
	}

	/**
	 * Refine the list of allowed mappings based on the current (partial)
	 * path (the start of the allowed mappings is set to the current solution and
	 * we prune the remainder, removing inconsistencies)
	 *
	 * @param allowedMappings
	 */
	void updateAllowedMappings(ListOfSets & allowedMappings) {
		bool mappingChanged = true;
		while (mappingChanged) {
			mappingChanged = false;
			size_t querySize = info.getSizeOfQuery();

			/*
			 * Look through each query node
			 * For each allowed queryNode that might map to a targetNode
			 * look at all the neighbors (queryAdjacentNode) to the query node and
			 * make sure that there is also an equivalent allowed neighbor in the
			 * target (targetAdjacentNode)
			 */
			for (size_t queryNode = 0; queryNode < querySize;
					queryNode++) {
				// need a copy here as allowedMappings.at(queryNode) is modified in loop
				auto queryNodes = allowedMappings.at(queryNode);
				for (size_t targetNode : queryNodes) {
					for (size_t queryAdjacentNode : queryAdjacencies.at(
							queryNode)) {
						bool match = false;
						for (size_t targetAdjacentNode : targetAdjacencies.at(
								targetNode)) {
							// if appropriate match edge type as well as just checking adjacency
							if (allowedMappings.at(queryAdjacentNode).count(
									targetAdjacentNode)
									&& info.sameEdgeType(queryNode,
											queryAdjacentNode, targetNode,
											targetAdjacentNode)) {
								match = true;
							}
						}
						if (!match) {
							allowedMappings.at(queryNode).erase(targetNode);
							mappingChanged = true;
						}
					}
				}
			}
		}
	}

	/**
	 * Builds query/target adjacency matrices
	 *
	 * @param query
	 */
	void createAdjacencies(bool query) {
		ListOfSets & matrix = query ? queryAdjacencies : targetAdjacencies;
		size_t size = query ? querySize : targetSize;
		for (size_t i = 0; i < size; i++) {
			for (size_t j = i + 1; j < size; j++) {
				bool adjacent =
						query ? info.isAdjacentInQuery(i, j) : info.isAdjacentInTarget(
										i, j);
				if (adjacent) {
					matrix.at(i).insert(j);
					matrix.at(j).insert(i);
				}
			}
		}
	}

};

} /* namespace GarethUtil */

#endif /* SRC_UTIL_ULLMAN_H_ */
