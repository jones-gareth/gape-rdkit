/*
 * MultiConformerGyrationShapeEvaluation.h
 *
 *  Created on: Apr 12, 2016
 *      Author: gjones
 */

#ifndef SRC_MOL_MULTICONFORMERGYRATIONSHAPEEVALUATION_H_
#define SRC_MOL_MULTICONFORMERGYRATIONSHAPEEVALUATION_H_

#include <memory>
#include "../util/Array2D.h"
#include "MulticonformerMolecule.h"

namespace GarethMol {

using namespace std;
using namespace GarethUtil;

class MultiConformerGyrationShapeEvalution {
public:
    MultiConformerGyrationShapeEvalution(const MulticonformerMolecule & mol):
        molecule(mol), distanceMatrix(mol.nConformers(), mol.nConformers())
    {
        buildDescriptorValues();
    }

    virtual ~MultiConformerGyrationShapeEvalution() {
    }

    MultiConformerGyrationShapeEvalution(const MultiConformerGyrationShapeEvalution & rhs) = delete;
    MultiConformerGyrationShapeEvalution & operator =(const MultiConformerGyrationShapeEvalution & rhs) = delete;
    MultiConformerGyrationShapeEvalution(MultiConformerGyrationShapeEvalution && rhs) = delete;
    MultiConformerGyrationShapeEvalution & operator =(MultiConformerGyrationShapeEvalution && rhs) = delete;

    const Array2D<double>& getDistanceMatrix() const {
        return distanceMatrix;
    }

    const int nItems() const {
        return molecule.nConformers();
    }

    const double itemDistance(const size_t i, const size_t j) const {
        return distanceMatrix.get(i, j);
    }

    const std::vector<std::vector<size_t> >  determineNearestNeighbourLists();

    const std::vector<std::vector<size_t> >& getNearestNeighbours() const {
        return nearestNeighbours;
    }

private:
    static const int nDescriptors = 4;
    const MulticonformerMolecule & molecule;
    Array2D<double> distanceMatrix;
    unique_ptr<Array2D<double>> descriptorValues = nullptr;
    std::vector<std::vector<size_t>> nearestNeighbours;

    void buildDescriptorValues();

};

} /* namespace GarethMol */

#endif /* SRC_MOL_MULTICONFORMERGYRATIONSHAPEEVALUATION_H_ */
