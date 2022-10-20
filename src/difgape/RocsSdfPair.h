/*
 * RocsSdfPair.h
 *
 *  Created on: May 16, 2014
 *      Author: Gareth Jones
 */

#ifndef ROCSSDFPAIR_H_
#define ROCSSDFPAIR_H_

#include <memory>
#include <string>

#include "../mol/mol.h"
#include "../mol/Molecule.h"
#include "../util/Array2D.h"
#include "RocsSdfMolecule.h"

namespace Difgape {

using namespace std;
using namespace GarethUtil;

class RocsSdfPair;
using RocsSdfPairPtr = shared_ptr<RocsSdfPair>;

class RocsSdfPair {

public:


	RocsSdfPair(RocsSdfMolecule::RocsSdfMoleculePtr q,
			RocsSdfMolecule::RocsSdfMoleculePtr r, string f);

	virtual ~RocsSdfPair() {
	}
	;

	double getDistance(int i, int j) {
		return distances(i, j);
	}

	string info() {
		return "Pair query " + queryName + " target " + targetName + " in file "
				+ file;
	}

	const double getScore(ScoreName scoreName) {
		return targetReference->getScore(scoreName);
	}


	const RocsSdfMolecule::RocsSdfMoleculePtr& getQueryReference() const {
		return queryReference;
	}

	const RocsSdfMolecule::RocsSdfMoleculePtr& getTargetReference() const {
		return targetReference;
	}

private:

	RocsSdfMolecule::RocsSdfMoleculePtr queryReference, targetReference;
	Array2D<double> distances;
	const string file, queryName, targetName;
	RocsSdfPair(const RocsSdfPair &) = delete;
	RocsSdfPair & operator =(const RocsSdfPair &) = delete;
	Molecule::MoleculePtr getMolecule(string name) const;

};

} /* namespace Difgape */

#endif /* ROCSSDFPAIR_H_ */
