/*
 * MultiConformerRmsEvalution.h
 *
 *  Created on: Sep 15, 2015
 *      Author: gjones
 */

#ifndef SRC_MOL_MULTICONFORMERRMSEVALUTION_H_
#define SRC_MOL_MULTICONFORMERRMSEVALUTION_H_

#include <stddef.h>
#include <vector>

#include "../util/Array2D.h"
#include "MulticonformerMolecule.h"

namespace GarethMol {

using namespace GarethUtil;
using namespace std;

/**
 * A class to compare all conformers in a multiconformer molecule and create an
 * rmsd matrix between conformers
 *
 * Rmsd are determined from the minimum value amongst all isomorphisms.  The isomorphisms
 * are only determined once.
 */
class MultiConformerRmsEvalution {
public:
	MultiConformerRmsEvalution(const MulticonformerMolecule & mol) :
			molecule(mol), rmsValues(mol.nConformers(), mol.nConformers()) {
	}

	virtual ~MultiConformerRmsEvalution() {
	}

	/**
	 * Determine all rms values and return array.
	 *
	 * @return
	 */
	const Array2D<double> & evaluateRms(bool matchTypes = true);

	/**
	 * Create nearest neighbour lists for all conformers
	 *
	 * @return
	 */
	const std::vector<std::vector<size_t>> & determineNearestNeighbourLists();

	/**
	 * Check the lists are correctly formed
	 *
	 * @return
	 */
	bool checkNeighbourLists() const;



	MultiConformerRmsEvalution(const MultiConformerRmsEvalution & rhs) = delete;
	MultiConformerRmsEvalution & operator =(
			const MultiConformerRmsEvalution & rhs) = delete;
	MultiConformerRmsEvalution(MultiConformerRmsEvalution && rhs) = delete;
	MultiConformerRmsEvalution & operator =(MultiConformerRmsEvalution && rhs) = delete;

	const std::vector<std::vector<size_t> >& getNearestNeighbours() const {
		return nearestNeighbours;
	}

	const Array2D<double>& getRmsValues() const {
		return rmsValues;
	}

	const size_t nItems() const {
	    return molecule.nConformers();
	}

	const double itemDistance(size_t i, size_t j) const {
	    return rmsValues.get(i, j);
	}
private:
	const MulticonformerMolecule & molecule;
	Array2D<double> rmsValues;
	// array of stored isomorphisms
	std::vector<std::vector<size_t>> isomorphisms;
	// list of nearest neighbours
	std::vector<std::vector<size_t>> nearestNeighbours;

	/**
	 * The callback function for the Ullman algorithm
	 * @param queryIdsToTargetIds
	 */
	void addIsomorphism(const std::vector<size_t> & queryIdsToTargetIds);


};

} /* namespace GarethUtil */

#endif /* SRC_MOL_MULTICONFORMERRMSEVALUTION_H_ */
