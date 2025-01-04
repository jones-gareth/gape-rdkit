#include <cassert>
#include <algorithm>
#include <set>
#include <boost/optional.hpp>
#include "Sssr.h"

namespace Gape {

    using namespace std;
    using namespace Gape;

/**
 * Check Paths to see if intermediate paths give a new or alternative shortest path.
 *
 * @param graphPath1
 * @param graphPath2
 */
    void GraphPath::checkShortestPath(const GraphPath &graphPath1,
                                      const GraphPath &graphPath2) {
        assert(graphPath1.getNode1() == node1);
        assert(graphPath2.getNode2() == node2);
        assert(graphPath1.getNode2() == graphPath2.getNode1());

        int testPathLength = graphPath1.getPathLength()
                             + graphPath2.getPathLength();
        if (testPathLength < pathLength) {
            REPORT(Reporter::TRACE) << "Found new shortest path between nodes "
                                    << node1 << " and " << node2 << " of length " << testPathLength;
            shortestPaths.clear();
            pathLength = testPathLength;
            addPaths(graphPath1, graphPath2, false);
        } else if (testPathLength == pathLength) {
            // alternative shortest path add to list
            REPORT(Reporter::TRACE)
                << "Found possible alternate shortest path between nodes "
                << node1 << " and " << node2 << " of length " << testPathLength;
            addPaths(graphPath1, graphPath2, false);
        }
    }

/**
 * Check to see the intermediate paths give a new or alternative path of one longer
 * than the shortest path between two nodes.
 *
 * @param graphPath1
 * @param graphPath2
 */
    void GraphPath::checkShortestPathPlus1(const GraphPath &graphPath1,
                                           const GraphPath &graphPath2) {
        assert(graphPath1.getNode1() == node1);
        assert(graphPath2.getNode2() == node2);

        int testPathLength = graphPath1.getPathLength()
                             + graphPath2.getPathLength();
        if (testPathLength == pathLength + 1) {
            REPORT(Reporter::TRACE)
                << "Found possible alternate shortest path plus one between nodes "
                << node1 << " and " << node2 << " of length " << testPathLength;
            addPaths(graphPath1, graphPath2, true);
        }

    }

/**
 * Join two intermediate paths
 *
 * @param path1
 * @param path2
 * @return
 */
    vector<int> GraphPath::joinPaths(const vector<int> &path1,
                                     const vector<int> &path2) {
        assert(*path1.cbegin() == node1);
        assert(*(path2.cend() - 1) == node2);
        assert(*(path1.cend() - 1) == *path2.cbegin());
        auto newPath = vector<int>(path1);
        newPath.reserve(path1.size() + path2.size());
        newPath.insert(newPath.cend(), path2.begin() + 1, path2.end());
        assert(newPath.size() == path1.size() + path2.size() - 1);
        assert(*newPath.cbegin() == node1);
        assert(*(newPath.end() - 1) == node2);
        return newPath;

    }

/**
 * Determine if a new path is unique, or if any of it's nodes are present in the current
 * path set (excluding start and end nodes)
 *
 * @param newPath
 * @param paths
 * @return
 */
    bool uniquePath(vector<int> newPath, vector<vector<int>> paths) {
        if (paths.size() == 0)
            return true;

        for (const auto &currentPath: paths) {
            for (auto iter = newPath.cbegin() + 1; iter < newPath.end() - 1;
                 ++iter) {
                if (find(currentPath.cbegin(), currentPath.cend(), *iter)
                    != currentPath.end()) {
                    return false;
                }
            }
        }
        return true;
    }

/**
 * Merge all possible intermediate paths.
 *
 * @param graphPath1
 * @param graphPath2
 * @param plus1
 */
    void GraphPath::addPaths(const GraphPath &graphPath1,
                             const GraphPath &graphPath2, const bool plus1) {
        assert(graphPath1.getNode1() == node1);
        assert(graphPath2.getNode2() == node2);

        for (const auto &path1: graphPath1.getShortestPaths()) {
            for (const auto &path2: graphPath2.getShortestPaths()) {
                const auto newPath = joinPaths(path1, path2);
                bool unique = true;
                if (!uniquePath(newPath, shortestPaths)) {
                    unique = false;
                }
                if (unique && plus1 && !uniquePath(newPath, shortestPathsPlus1)) {
                    unique = false;
                }
                if (unique) {
                    if (plus1) {
                        shortestPathsPlus1.push_back(newPath);
                    } else {
                        shortestPaths.push_back(newPath);
                    }
                }
            }
        }
    }

/**
 * Creates ring from two alternate paths between nodes
 *
 * @param segment1
 * @param segment2
 * @return
 */
    vector<int> joinRingSegments(const vector<int> &segment1,
                                 const vector<int> &segment2) {
        assert(*segment1.cbegin() == *segment2.cbegin());
        assert(*(segment1.cend() - 1) == *(segment2.cend() - 1));
        vector<int> segment2Reverse(segment2.size());
        //segment2Reverse.reserve(segment2.size());
        reverse_copy(segment2.begin(), segment2.end(), segment2Reverse.begin());
        auto newPath = vector<int>(segment1);
        newPath.reserve(segment1.size() + segment2.size());
        //newPath.insert(newPath.cbegin(), segment1.cbegin(), segment1.cend());
        newPath.insert(newPath.cend(), segment2Reverse.cbegin() + 1,
                       segment2Reverse.cend() - 1);
        return newPath;
    }

/**
 * Return all possible rings between two nodes
 *
 * @return
 */
    vector<vector<int>> GraphPath::toRings() const {
        if (pathLength < 1) {
            return {};
        }
        vector<vector<int>> possibleRings = {};
        for (size_t i = 0; i < shortestPaths.size(); i++) {
            for (size_t j = i + 1; j < shortestPaths.size(); j++) {
                auto possibleRing = joinRingSegments(shortestPaths.at(i),
                                                     shortestPaths.at(j));
                possibleRings.push_back(possibleRing);
            }
        }
        for (size_t i = 0; i < shortestPaths.size(); i++) {
            for (size_t j = 0; j < shortestPathsPlus1.size(); j++) {
                auto possibleRing = joinRingSegments(shortestPaths.at(i),
                                                     shortestPathsPlus1.at(j));
                possibleRings.push_back(possibleRing);
            }
        }
        return possibleRings;
    }

/**
 * Test to see if a ring has already been identified, and add it to the set if it hasn't
 *
 * @param newRing
 * @param rings
 * @return
 */
    bool checkRing(const vector<int> &newRing, vector<vector<int>> &rings) {
        REPORT(Reporter::TRACE) << "Testing ring " << collectionToString(newRing)
                                << " for presence";
        for (const auto &ring: rings) {
            if (newRing.size() == ring.size()) {
                bool sameRing = true;
                for (const auto node: newRing) {
                    if (find(ring.cbegin(), ring.cend(), node) == ring.cend()) {
                        sameRing = false;
                        break;
                    }
                }
                if (sameRing) {
                    REPORT(Reporter::TRACE) << "Ring already present";
                    return false;
                }
            }
        }
        REPORT(Reporter::TRACE) << "New ring found";
        rings.push_back(newRing);
        return true;
    }

/**
 *  Builds an SSSR from all unique rings.
 *
 *  The orgininal definition of the first nSssrRings (= nEdges - nNodes + nFragments)
 *  rings from the sorted by size list of rings is too simplistic. Consider a cube or tetrahedron.
 *
 *  This sorts the rings in order and rings are added to the growing SSSR, provided they are
 *  not completely covered by SSSR rings of a smaller size.
 *
 * @param allRings
 * @return
 */
    const vector<vector<int>> extractSssr(vector<vector<int>> &allRings) {

        // sort rings by size
        auto sortFn =
                [](const vector<int> &a, const vector<int> &b) { return a.size() < b.size(); };
        sort(allRings.begin(), allRings.end(), sortFn);
        vector<vector<int>> sssrRings = {};

        set<int> coveredNodes = {};
        size_t currentRingSize = allRings.at(0).size();

        for (auto ring: allRings) {

            if (ring.size() > currentRingSize) {
                // moved to a new size - update set of covered nodes
                for (auto ring2: sssrRings) {
                    if (ring2.size() > currentRingSize) {
                        break;
                    }
                    for (int node: ring2) {
                        coveredNodes.insert(node);
                    }
                }
                currentRingSize = ring.size();
            }

            // see if ring is covered by smaller rings
            bool covered = true;
            for (int node: ring) {
                if (coveredNodes.count(node) == 0) {
                    covered = false;
                    break;
                }
            }

            if (!covered) {
                sssrRings.push_back(ring);
            }

        }

        return sssrRings;

    }

} /* namespace Gape */
