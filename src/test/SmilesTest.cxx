/*
 * loadMol.cpp
 *
 *	Test file for reading in smiles
 *
 *  Created on: Apr 29, 2014
 *      Author: gjones
 */

#include <gtest/gtest.h>
#include <stddef.h>
#include <cstdio>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "../mol/Molecule.h"
#include "../mol/MolComparer.h"
#include "../mol/SmilesParser.h"
#include "../util/Reporter.h"
#include "../util/FileSystemUtil.h"

namespace {

using namespace std;
using namespace GarethUtil;
using namespace GarethMol;
using namespace ::testing;

/**
 * Test class for reading molecules
 */
class SmilesTest: public Test {

protected:
	SmilesTest() {
		Reporter::setMinReportingLevel(Reporter::NORMAL);
	}
};

bool sameMolecule(const Molecule & moleculeA, const Molecule & moleculeB) {
	MolComparer molComparer(moleculeA, moleculeB);
	molComparer.setSubgraph(false);
	molComparer.setHeavyAtomOnly(false);
	auto callback = [](const vector<size_t> & queryIdsToTargetIds) {;};
	molComparer.setCallbackFunction(callback);
	molComparer.compare();
	return molComparer.getIsomorphisms() > 0;
}

TEST_F(SmilesTest, testEquiv1) {
	SmilesParser parser;
	auto smiles1 = "CC1=CC(Br)CCC1";
	auto smiles2 = "CC1=CC(CCC1)Br";
	auto moleculeA = parser.parseSmiles(smiles1);
	auto moleculeB = parser.parseSmiles(smiles2);
	auto equiv = sameMolecule(*moleculeA, *moleculeB);
	ASSERT_TRUE(equiv);

	parser.parseSmiles("C(=O)Oc1c(F)c(F)c(F)c(F)c1(F)");
}

TEST_F(SmilesTest, testError) {
	SmilesParser parser;
	parser.parseSmiles("CCCCSCC 19");
}

bool readSmilesFile(istream & in) {
	SmilesParser parser;
	auto count = 0;
	while (!in.eof()) {
		string line;
		getline(in, line);
		removeTrailingLF(line);
		trim(line);
		if (line == "") {
			continue;
		}
		if (startsWith(line, "#")) {
			continue;
		}
		auto smiles = line;
		REPORT(Reporter::DETAIL) << "Parsing smiles " << smiles;
		auto molecule = parser.parseSmiles(smiles);
		molecule->initialize();
		count++;
		if (count % 20 == 0) {
			REPORT(Reporter::INFO) << "Parsed "<<count<<" smiles";
		}
	}
	return true;
}

TEST_F(SmilesTest,  testNci) {
	//string fileName("../resources/NCI-Open_09-03.smi.gz");
	string fileName("../resources/NCI-Open_09-03_20K.smi");
	readObjectFromFile<bool>(fileName, readSmilesFile);
}


}
