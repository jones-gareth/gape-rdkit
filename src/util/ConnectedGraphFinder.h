/*
 * ConnectedGraphFinder.h
 *
 *  Created on: Nov 5, 2015
 *      Author: gjones
 */

#ifndef SRC_UTIL_CONNECTEDGRAPHFINDER_H_
#define SRC_UTIL_CONNECTEDGRAPHFINDER_H_

namespace Gape {

#include <vector>
#include <cassert>

using namespace std;

/**
 * A class for determine fully connected sub-graphs within a larger graph
 *
 * Instantiate with a policy class that has two methods
 *
 * Returns the total number of nodes in the graph
 * const int getNumberNodes() const;
 *
 * Returns nodes connected to a given node
 * const vector<int> getNodeNeighbours(const int nodeNo) const;
 *
 */
template<typename GraphInfo>
class ConnectedGraphFinder {
public:
	ConnectedGraphFinder(const GraphInfo & graphInfo_) :
			graphInfo(graphInfo_), fragmentNos(graphInfo.getNumberNodes(), -1) {
	}

	virtual ~ConnectedGraphFinder() {
	}

	ConnectedGraphFinder(const ConnectedGraphFinder & rhs) = delete;
	ConnectedGraphFinder & operator =(const ConnectedGraphFinder & rhs) = delete;
	ConnectedGraphFinder(ConnectedGraphFinder && rhs) = delete;
	ConnectedGraphFinder & operator =(ConnectedGraphFinder && rhs) = delete;

	const vector<int> & findConnectedGraphs() {
		nFragments = 0;
		for (int nodeNo = 0; nodeNo < graphInfo.getNumberNodes(); nodeNo++) {
			if (fragmentNos.at(nodeNo) == -1) {
				findFragment(nodeNo, nFragments);
				nFragments++;
			}
		}

		return fragmentNos;
	}

	/**
	 * Returns a vector with position being node number and value being graph number.
	 *
	 * @return
	 */
	const vector<int>& getFragmentNos() const {
		return fragmentNos;
	}

	/**
	 * Returns the number of fully connected graphs
	 * @return
	 */
	int getNFragments() const {
		return nFragments;
	}

private:
	const GraphInfo & graphInfo;
	vector<int> fragmentNos;
	int nFragments = 0;

	void findFragment(int nodeNo, int fragmentNo) {
		assert(fragmentNos.at(nodeNo) == -1);
		fragmentNos.at(nodeNo) = fragmentNo;

		for (auto neighbourNodeNo : graphInfo.getNodeNeighbours(nodeNo)) {
			int neighbourFragmentNo = fragmentNos.at(neighbourNodeNo);
			if (neighbourFragmentNo == -1) {
				findFragment(neighbourNodeNo, fragmentNo);
			} else {
				assert(neighbourFragmentNo == fragmentNo);
			}
		}
	}

};

} /* namespace Gape */

#endif /* SRC_UTIL_CONNECTEDGRAPHFINDER_H_ */
