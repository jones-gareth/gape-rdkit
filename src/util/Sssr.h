/*
 * Sssr.h
 *
 *  Created on: Oct 30, 2015
 *      Author: Gareth Jones
 */

#ifndef SRC_UTIL_SSSR_H_
#define SRC_UTIL_SSSR_H_

#include "Array2D.h"
#include "Reporter.h"
#include "ConnectedGraphFinder.h"
#include <vector>

namespace Gape {

using namespace std;

/**
 * Class to represent paths between two nodes in a graph.
 *
 * Stores unique shortest paths (alternative paths cannot intersect) and shortest
 * paths plus one (they must also be unique and cannot interest with any of the
 * shortest paths). The intersect criteria is relevent to detection of rings.
 */
class GraphPath {

public:
	GraphPath(int node1_, int node2_, bool edge, int max) :
			node1(node1_), node2(node2_) {
		shortestPaths = {};
		shortestPathsPlus1 = {};
		// initialize paths
		if (node1 == node2) {
			pathLength = 0;
		} else if (edge) {
			REPORT(Reporter::TRACE) << "Edge between nodes "<<node1<<" and "<<node2;
			pathLength = 1;
			shortestPaths = { {node1, node2}};
		} else {
			pathLength = max;
		}
	}

	virtual ~GraphPath() {
	}

	GraphPath(const GraphPath & rhs) = delete;
	GraphPath & operator =(const GraphPath & rhs) = delete;
	GraphPath(GraphPath && rhs) = delete;
	GraphPath & operator =(GraphPath && rhs) = delete;

	const int getNode1() const {
		return node1;
	}

	const int getNode2() const {
		return node2;
	}

	int getPathLength() const {
		return pathLength;
	}

	const vector<vector<int> >& getShortestPaths() const
	{
		return shortestPaths;
	}

	const vector<vector<int> >& getShortestPathsPlus1() const
	{
		return shortestPathsPlus1;
	}

	/**
	 * Check Paths to see if intermediate paths give a new or alternative shortest path.
	 *
	 * @param graphPath1
	 * @param graphPath2
	 */
	void checkShortestPath(const GraphPath & graphPath1, const GraphPath & graphPath2);

	/**
	 * Check to see the intermediate paths give a new or alternative path of one longer
	 * than the shortest path between two nodes.
	 *
	 * @param graphPath1
	 * @param graphPath2
	 */
	void checkShortestPathPlus1(const GraphPath & graphPath1, const GraphPath & graphPath2);

	/**
	 * See if the path information between nodes describes a ring
	 *
	 * @return
	 */
	vector<vector<int>> toRings() const;

private:
	const int node1, node2;
	vector<vector<int>> shortestPaths, shortestPathsPlus1;
	int pathLength=0;

	/**
	 * Join two intermediate paths
	 *
	 * @param path1
	 * @param path2
	 * @return
	 */
	vector<int> joinPaths(const vector<int> & path1, const vector<int> & path2);

	/**
	 * Merge all possible intermediate paths.
	 *
	 * @param graphPath1
	 * @param graphPath2
	 * @param plus1
	 */
	void addPaths(const GraphPath & graphPath1, const GraphPath & graphPath2, const bool plus1);

};

bool checkRing(const vector<int> & ring, vector<vector<int>> & rings);
const vector<vector<int>> extractSssr(vector<vector<int>> & allRings);

/**
 * Determines SSSR for a graph.  Implements the method described in
 *
 * 1.Lee, C. J., Kang, Y.-M., Cho, K.-H. & No, K. T.
 * A robust method for searching the smallest set of smallest rings with a path-included distance matrix.
 * PNAS 106, 17355–17358 (2009).
 *
 * Should be instantiated with a policy class with the following methods:
 *
 * Returns the total number of nodes in the graph
 * const int getNumberNodes() const;
 *
 * Returns nodes connected to a given node
 * const vector<int> getNodeNeighbours(const int nodeNo) const;
 *
 * Returns true if the two nodes are connected by an edge
 * bool isEdge(const int & node1, const int & node2) const
 */
template<typename GraphInfo>
class Sssr {
public:
	Sssr(const GraphInfo & graphInfo_) :
			graphInfo(graphInfo_), graphPaths(graphInfo.getNumberNodes(),
					graphInfo.getNumberNodes()) {
	}

	virtual ~Sssr() {
	}

	Sssr(const Sssr & rhs) = delete;
	Sssr & operator =(const Sssr & rhs) = delete;
	Sssr(Sssr && rhs) = delete;
	Sssr & operator =(Sssr && rhs) = delete;

	/**
	 * Determine the SSSR.
	 *
	 * @return
	 */
	const vector<vector<int>> & doSssr() {
		const int nNodes = graphInfo.getNumberNodes();
		int nEdges = 0;

		// initialize all path information
		for (auto node1 = 0; node1 < nNodes; node1++) {
			for (auto node2 = 0; node2 < nNodes; node2++) {
				const bool isEdge = graphInfo.isEdge(node1, node2);
				if (isEdge) {
					nEdges++;
				}
				auto graphPath = make_unique<GraphPath>(node1, node2, isEdge,
						nNodes);
				graphPaths(node1, node2) = move(graphPath);
			}
		}
		nEdges = nEdges / 2;

		// Determine size of SSSR
		ConnectedGraphFinder<GraphInfo> connectedGraphFinder(graphInfo);
		connectedGraphFinder.findConnectedGraphs();
		int nFragments = connectedGraphFinder.getNFragments();
		int nSssrRings = nEdges - nNodes + nFragments;
		REPORT(Reporter::TRACE) << "# edges " << nEdges << " # nodes "
				<< nNodes << " # sssr " << nSssrRings;
		if (nSssrRings == 0) {
			rings = {};
			return rings;
		}

		// Determine shortest paths using Floyd–Warshall algorithm
		for (auto k = 0; k < nNodes; k++) {
			for (auto node1 = 0; node1 < nNodes; node1++) {
				for (auto node2 = 0; node2 < nNodes; node2++) {
					auto & graphPath = graphPaths(node1, node2);
					auto &path1 = graphPaths(node1, k);
					auto &path2 = graphPaths(k, node2);
					graphPath->checkShortestPath(*path1, *path2);
				}
			}
		}

		// Determine paths one longer than shortest paths
		for (auto k = 0; k < nNodes; k++) {
			for (auto node1 = 0; node1 < nNodes; node1++) {
				for (auto node2 = 0; node2 < nNodes; node2++) {
					auto & graphPath = graphPaths(node1, node2);
					auto &path1 = graphPaths(node1, k);
					auto &path2 = graphPaths(k, node2);
					graphPath->checkShortestPathPlus1(*path1, *path2);
				}
			}
		}

		// Extract rings from path information
		vector<vector<int>> allRings = { };
		for (auto node1 = 0; node1 < nNodes; node1++) {
			for (auto node2 = 0; node2 < nNodes; node2++) {
				auto & graphPath = graphPaths(node1, node2);
				auto possibleRings = graphPath->toRings();
				for (auto ring : possibleRings) {
					if (checkRing(ring, allRings)) {
						REPORT(Reporter::TRACE) << "Found unique ring "
								<< collectionToString(ring);
					}
				}
			}
		}

		if (false) {
			// Extract SSSR - this original definition is too simplistic- see this smile and consider
			// cubes and tetrahedrons.
			// c1cccc(c1)C23#C4([Co]2([Co]34([CH-]#[O+])([CH-]#[O+])[CH-]#[O+])([CH-]#[O+])([CH-]#[O+])[CH-]#[O+])c5ccccc5

			// Sort all rings by size
			assert(static_cast<int>(allRings.size()) >= nSssrRings);
			auto sortFn =
					[] (const vector<int> &a, const vector<int> &b) {return a.size() < b.size();};
			sort(allRings.begin(), allRings.end(), sortFn);
			rings = {};

			rings.insert(rings.begin(), allRings.begin(),
					allRings.begin() + nSssrRings);
		} else {
			rings = extractSssr(allRings);
			if (rings.size() != static_cast<size_t>(nSssrRings)) {
				REPORT(Reporter::DETAIL) << "SSSR Expected " <<nSssrRings << " rings but found " << rings.size();
			}
		}

		REPORT(Reporter::TRACE) << "Sssr determined:";
		for (auto ring : rings) {
			REPORT(Reporter::TRACE) << collectionToString(ring);
		}

		return rings;
	}

	const vector<vector<int> >& getRings() const {
		return rings;
	}

private:
	const GraphInfo & graphInfo;
	//UniquePtrArray2D<GraphPath> graphPaths;
	Array2D<unique_ptr<GraphPath>> graphPaths;
	vector<vector<int>> rings = { };

};

} /* namespace Gape */

#endif /* SRC_UTIL_SSSR_H_ */
