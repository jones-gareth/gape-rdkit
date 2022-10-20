/*
 * loadMol.cpp
 *
 *	Test file for reading in molecules
 *
 *  Created on: Apr 29, 2014
 *      Author: gjones
 */


#include <cassert>
#include <iostream>
#include <string>
#include <vector>

#include "../mol/Molecule.h"
#include "../util/Reporter.h"


using namespace std;
using namespace GarethUtil;
using namespace GarethMol;

/**
 * Test routine for reading molecules
 *
 * @param argc
 * @param argv
 * @return
 */
int main(int argc, char* argv[]) {
	assert(argc == 3);
	string inFile = string(argv[1]);
	string outFile = string(argv[2]);

	Reporter::setMinReportingLevel(Reporter::NORMAL);
	REPORT(Reporter::INFO) << "loading file " << inFile;

	vector<Molecule::MoleculePtr> mols = Molecule::readMoleculesFromFile(inFile);
	Molecule::writeMoleculesToFile(outFile, mols);

	return 1;
}

