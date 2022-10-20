/*
 * MolecularRms.h
 *
 *  Created on: Sep 10, 2015
 *      Author: gjones
 */

#ifndef SRC_MOL_MOLECULARRMS_H_
#define SRC_MOL_MOLECULARRMS_H_

#include "Molecule.h"
#include "MolComparer.h"

namespace GarethMol {

using namespace std;
using namespace GarethUtil;

/**
 * A class for comparing two structures and creating a smart rms
 */
class MolecularRms {
public:
	MolecularRms(const Molecule & molA_, const Molecule & molB_) :
			queryMolecule(molA_), targetMolecule(molB_) {

	}

	/**
	 * Calculate the rms (accounting for isomorphisms) between the coordinates of
	 * the constructor molecules.  Calls determineIsomorphisms() prior to calculating RMS
	 *
	 * @return
	 */
	double determineRms() {
		determineIsomorphisms();
		return determineRms(queryMolecule.getCoords(),
				targetMolecule.getCoords());
	}

	/**
	 * Calculate the rms (accounting for isomorphisms between the constructor molecules),
	 * but using the suppled coordinates.  Must call determineIsomorphisms() prior to
	 * using this routine.  For use in analysing multiconformer molecules.
	 *
	 * @param queryCoords
	 * @param targetCoords
	 * @return
	 */
	double determineRms(const CoordMatrix & queryCoords,
			const CoordMatrix & targetCoords);

	/**
	 * Finds all isomorphisms between the query and target molecule
	 */
	void determineIsomorphisms();

	virtual ~MolecularRms();

	/**
	 * Create a comparer that compares a molecule with itself
	 *
	 * @param molecule
	 * @return
	 */
	static unique_ptr<MolecularRms> createSelfComparer(
			const Molecule &molecule);

	MolecularRms(const MolecularRms & rhs) = delete;
	MolecularRms & operator =(const MolecularRms & rhs) = delete;
	MolecularRms(MolecularRms && rhs) = delete;
	MolecularRms & operator =(MolecularRms && rhs) = delete;

	/**
	 * Return all the isomorphisms found
	 * @return
	 */
	const std::vector<double>& getIsomorphismRmsds() const {
		return isomorphismRmsds;
	}

	/**
	 * Return the minimum rmsd found.
	 * @return
	 */
	double getMinRmsd() const {
		return minRmsd;
	}

	/**
	 * A static method that does the least squares fitting between two sets of atomic coordinates.
	 *
	 * @param queryIdsToTargetIds
	 * @param queryAtoms
	 * @param targetAtoms
	 * @param queryCoords
	 * @param targetCoords
	 * @param doLeastSquaresFit
	 * @return
	 */
	static double leastSquaresAtomsFit(
			const std::vector<size_t> & queryIdsToTargetIds,
			const std::vector<const Atom *> & queryAtoms,
			const std::vector<const Atom *> & targetAtoms,
			const CoordMatrix & queryCoords, const CoordMatrix & targetCoords,
			bool doLeastSquaresFit = true);

	bool isDoLeastSquaresFit() const {
		return doLeastSquaresFit;
	}

	void setDoLeastSquaresFit(bool doLeastSquaresFit) {
		this->doLeastSquaresFit = doLeastSquaresFit;
	}

private:
	const Molecule & queryMolecule, &targetMolecule;
	bool doLeastSquaresFit { true };
	CoordMatrix targetCoords, queryCoords;
	unique_ptr<MolComparer> molComparer;
	void addIsomorphism(const std::vector<size_t> & queryIdsToTargetIds);
	void evaluateIsomorphisms();
	std::vector<std::vector<size_t>> isomorphisms;

	double minRmsd = numeric_limits<double>::max();
	std::vector<double> isomorphismRmsds { };
};

} /* namespace GarethUtil */

#endif /* SRC_MOL_MOLECULARRMS_H_ */
