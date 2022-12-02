/*
 * Array2D.cpp
 *
 *  Created on: Feb 20, 2016
 *      Author: gjones
 */

#include "Array2D.h"
#include "Reporter.h"

namespace Gape {

template<>
void Array2D<double>::normalizeColumns() {

    for (auto columnNo = 0ul; columnNo < nColumns; columnNo++) {
        double sum = .0;
        double sqrSum = .0;

        for (auto rowNo = 0ul; rowNo < nRows; rowNo++) {
            auto val = get(rowNo, columnNo);
            sum += val;
            sqrSum += val * val;
        }

        auto dSize = static_cast<double>(nRows);
        auto mean = sum / dSize;
        auto sumOfSquares = sqrSum - (sum * sum) / dSize;
        auto variance = sumOfSquares / (dSize - 1.0);
        auto stdDev = sqrt(variance);

        REPORT(Reporter::NORMAL) << "column " << columnNo << " mean " << mean
                << " std deviation " << stdDev;

        for (auto rowNo = 0ul; rowNo < nRows; rowNo++) {
            auto val = get(rowNo, columnNo);
            auto z = (val - mean) / stdDev;
            set(rowNo, columnNo, z);
        }

    }

}

template<>
double Array2D<double>::rowSqrDistance(const int i, const int j) {
    double sum = .0;
    for (auto c = 0ul; c < nColumns; c++) {
        auto diff = get(i, c) - get(j, c);
        sum += diff * diff;
    }
    return sum;
}

}
