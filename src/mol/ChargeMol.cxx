/*
 * ChargeMol.cxx
 *
 *  Created on: Jan 21, 2016
 *      Author: gjones
 */

#include "Molecule.h"
#include "Charge.h"
#include "../util/Reporter.h"


using namespace GarethMol;
using namespace GarethUtil;

int main(int argc, char** argv) {
	Reporter::setMinReportingLevel(Reporter::DEBUG);
	Charge charge;

	if (argc != 3) {
		cerr << "Usage: "<<argv[0] << " <inFile> <outFile>";
		exit(EXIT_FAILURE);
	}

	auto molecules = Molecule::readMoleculesFromFile(argv[1]);
	for (auto & molecule : molecules) {
		molecule->initialize();
		REPORT(Reporter::INFO) << "Adding formal charges to molecule " << molecule->getName();
		charge.chargeMolecule(*molecule);
	}
	Molecule::writeMoleculesToFile(argv[2], molecules);

}
