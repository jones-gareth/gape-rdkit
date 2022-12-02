/*
 * DetermineNearestNeighbourLists.h
 *
 *  Created on: Apr 21, 2016
 *      Author: gjones
 */

#ifndef SRC_UTIL_DETERMINENEARESTNEIGHBOURLISTS_H_
#define SRC_UTIL_DETERMINENEARESTNEIGHBOURLISTS_H_

#include <vector>
#include <numeric>

namespace Gape {

#include "Reporter.h"

/**
 * A template class to determine nearest neighbour lists.
 *
 * The policy class NeighbourInfo needs to provide the following class methods
 *
 * size_t nItems() const - number of items
 * double itemDistance(size_t i, size_t j) const - distance between items i and j
 */
template<typename NeighbourInfo>
class DetermineNearestNeighbourLists {
public:

    DetermineNearestNeighbourLists(const DetermineNearestNeighbourLists & rhs) = delete;
    DetermineNearestNeighbourLists & operator =(
            const DetermineNearestNeighbourLists & rhs) = delete;
    DetermineNearestNeighbourLists(DetermineNearestNeighbourLists && rhs) = delete;
    DetermineNearestNeighbourLists & operator =(
            DetermineNearestNeighbourLists && rhs) = delete;

    /**
     * Creates nearest neighbour lists of indexes.
     *
     * @param info
     * @return
     */
    static std::vector<std::vector<size_t>> determineNearestNeighbourLists(
            const NeighbourInfo & info) {
        auto nearestNeighbours = std::vector<std::vector<size_t>>(
                info.nItems());

        size_t itemNo = 0;
        for (auto iter = nearestNeighbours.begin();
                iter < nearestNeighbours.end(); ++iter, itemNo++) {
            iter->resize(info.nItems());
            std::iota(iter->begin(), iter->end(), 0);

            auto sortFn =
                    [itemNo, & info](size_t no1, size_t no2) {
                        return info.itemDistance(itemNo, no1) < info.itemDistance(itemNo, no2);
                    };
            sort(iter->begin(), iter->end(), sortFn);
        }

        assert(checkNeighbourLists(info, nearestNeighbours));
        return nearestNeighbours;
    }

    /**
     * Checks ordering of neighbour lists is correct
     *
     * @param info
     * @param nearestNeighbours
     * @return
     */
    static bool checkNeighbourLists(const NeighbourInfo & info,
           const std::vector<std::vector<size_t>> & nearestNeighbours) {
        size_t nItems = info.nItems();
        // loop over all conformer lists
        size_t itemNo = 0;
        if (nearestNeighbours.size() != nItems) {
            REPORT(Reporter::WARN) << "Top level neighbour list size error!";
            return false;
        }

        for (auto iter = nearestNeighbours.begin();
                iter < nearestNeighbours.end(); ++iter, itemNo++) {
            if (iter->size() != nItems) {
                REPORT(Reporter::WARN) << "Neighbour list size error!";
                return false;
            }
            // check each conformer list is ordered via rms
            for (auto listIter = iter->begin(); listIter < iter->end();
                    ++listIter) {
                auto nextIter = next(listIter);
                if (nextIter != iter->end()) {
                    assert(*listIter != *nextIter);
                    auto rms = info.itemDistance(itemNo, *listIter);
                    auto next = info.itemDistance(itemNo, *nextIter);
                    if (rms > next) {
                        REPORT(Reporter::WARN)
                                << "Neighbour list sequence check failed!";
                        return false;
                    }
                }
            }
        }

        REPORT(Reporter::DEBUG) << "Neighbour list sequence check passed";
        return true;
    }

private:
    DetermineNearestNeighbourLists() {
    }

    virtual ~DetermineNearestNeighbourLists() {
    }

};

} /* namespace Gape */

#endif /* SRC_UTIL_DETERMINENEARESTNEIGHBOURLISTS_H_ */
