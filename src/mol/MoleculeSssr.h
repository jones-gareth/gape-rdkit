/*
 * MoleculeSssr.h
 *
 *  Created on: Nov 3, 2015
 *      Author: gjones
 */

#ifndef SRC_MOL_MOLECULESSSR_H_
#define SRC_MOL_MOLECULESSSR_H_

#include "../util/Sssr.h"
#include "Molecule.h"

using namespace std;
using namespace GarethUtil;

namespace GarethMol {

/**
 * Class to determine Smallest Set of Smallest Rings for a molecule.  Atoms do not have
 * to have their neighbourhood determined already.
 *
 * See Sssr.h template class
 */
class MoleculeSssr {
	friend Sssr<MoleculeSssr> ;
	friend ConnectedGraphFinder<MoleculeSssr> ;

public:
	MoleculeSssr(Molecule & molecule_) :
			molecule(molecule_) {
	}

	virtual ~MoleculeSssr() {
	}

	MoleculeSssr(const MoleculeSssr & rhs) = delete;
	MoleculeSssr & operator =(const MoleculeSssr & rhs) = delete;
	MoleculeSssr(MoleculeSssr && rhs) = delete;
	MoleculeSssr & operator =(MoleculeSssr && rhs) = delete;

	/**
	 * Build the SSSR
	 */
	void doSssr();


private:
	Molecule &molecule;

	/**
	 * Return number of atoms in molecule
	 * @return
	 */
	const int getNumberNodes() const {
		return molecule.nAtoms();
	}

	/**
	 * Return the atom indices of neighbouring atoms
	 * @param atomNo
	 * @return
	 */
	const vector<int> getNodeNeighbours(const int atomNo) const;

	/**
	 * Return true if two atoms are connected
	 *
	 * @param atomNo1
	 * @param atomNo2
	 * @return
	 */
	bool isEdge(const int & atomNo1, const int & atomNo2) const {
		return molecule.getBond(atomNo1, atomNo2) != nullptr;
	}
};

} /* namespace GarethMol */

#endif /* SRC_MOL_MOLECULESSSR_H_ */
