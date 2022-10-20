/*
 * RocsSdfPair.cpp
 *
 *  Created on: May 16, 2014
 *      Author: Gareth Jones
 */

#include "RocsSdfPair.h"

#include <fstream>

#include "../mol/mol.h"
#include "../util/CoordOps.h"
#include "../util/FileSystemUtil.h"

namespace Difgape {

using namespace std;
using namespace GarethUtil;

static void createDistanceMatrix(RocsSdfMolecule::RocsSdfMoleculePtr q,
		RocsSdfMolecule::RocsSdfMoleculePtr t,
		Array2D<double> & distances) {

	int nQueryFeatures = q->getNFeatures();
	int nTargetFeatures = t->getNFeatures();
	for (int i = 0; i < nQueryFeatures; i++) {
		for (int j = 0; j < nTargetFeatures; j++) {
			CoordVector qf = q->getFeatureCoordinate(i);
			CoordVector tf = t->getFeatureCoordinate(j);
			distances(i, j) = distance(qf, tf);
		}
	}
}

RocsSdfPair::RocsSdfPair(RocsSdfMolecule::RocsSdfMoleculePtr q,
		RocsSdfMolecule::RocsSdfMoleculePtr t, string f) :
		queryReference { q }, targetReference { t },
		distances {q->getNFeatures(), t->getNFeatures()},
		file { f }, queryName {
				q->getName() }, targetName { t->getName() } {
	createDistanceMatrix(q, t, distances);
}

Molecule::MoleculePtr RocsSdfPair::getMolecule(string name) const {

}

} /* namespace Difgape */
