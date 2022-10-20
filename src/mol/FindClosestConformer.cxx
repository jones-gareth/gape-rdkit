/*
 *
 * Finds the closest conformer to a structure
 *
 * FindClosestConformer.cpp
 *
 *  Created on: Feb 9, 2016
 *      Author: gjones
 */

#include <iostream>
#include <string>

#include <boost/format.hpp>

#include "Molecule.h"
#include "MulticonformerMolecule.h"
#include "MolComparer.h"
#include "MolecularRms.h"

using namespace std;
using namespace GarethMol;

int main(int argc, char* argv[]) {

	auto usage = (boost::format(
			"usage: %s <single conformer file> <multi conformer file>\n")
			% argv[0]).str();

	int c;

	while ((c = getopt(argc, argv, "")) != -1) {
		switch (c) {
		case '?':
			cout << usage;
			return EXIT_SUCCESS;
		default:
			cout << "Unknown option" << endl;
			cout << usage;
			return EXIT_FAILURE;
		}
	}

	if (argc - optind != 2) {
		cout << usage;
		return EXIT_FAILURE;
	}

	string queryFile(argv[optind]);
	string targetFile(argv[optind + 1]);

	auto molecule = Molecule::readAndInitializeMolecule(queryFile);
	auto multiMolecule = MulticonformerMolecule::readAndInitializeMolecule(
			targetFile);

	const auto & coords1 = molecule->getCoords();
	auto bestRms = numeric_limits<double>::max();
	int bestConformerNo = -1;
	MolecularRms molecularRms(*molecule, *multiMolecule);
	molecularRms.determineIsomorphisms();

	for (const auto & conformer : multiMolecule->getConformers()) {

		auto rmsd = molecularRms.determineRms(molecule->getCoords(),
				conformer->getCoordinates());
		REPORT(Reporter::DEBUG) << "Rmsd between ligand and conformer "
				<< to_string(conformer->getConformerNo() + 1) << ": " << rmsd;
		if (rmsd < bestRms) {
			bestRms = rmsd;
			bestConformerNo = conformer->getConformerNo();
		}
	}

	cout << "N_conformers "<<multiMolecule->nConformers()<< " closest " << to_string(bestConformerNo + 1)
			<< " rms " << bestRms << endl;

	return EXIT_SUCCESS;
}

