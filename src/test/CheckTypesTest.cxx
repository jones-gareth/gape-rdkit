/*
 * CheckTypesTest.cxx
 *
 *	Test file for automatic assignment of atom types
 *
 *  Needs some more tests
 *
 */

#include <gtest/gtest.h>
#include <stddef.h>
#include <cstdio>
#include <iostream>
#include <memory>
#include <string>
#include <vector>
#include <fstream>

#include "../mol/Molecule.h"
#include "../mol/MoleculeSssr.h"

namespace {

using namespace std;
using namespace GarethUtil;
using namespace GarethMol;
using namespace ::testing;

/**
 * Test class for automatic assignment of atom types
 */
class CheckTypesTest: public Test {
protected:

	CheckTypesTest() {
		Reporter::setMinReportingLevel(Reporter::NORMAL);
	}
};

/**
 * Check P assignment bug
 */
TEST_F(CheckTypesTest, test1r1h) {
	auto ligandFile = "../resources/1r1h_ligand.mol";
	Molecule molecule;
	assert(molecule.readMoleculeFromFile(ligandFile));
	MoleculeSssr moleculeSssr(molecule);
	moleculeSssr.doSssr();
	molecule.createAtomNeighbourhood();
	int typesChanged = molecule.assignAtomTypes();
	REPORT(Reporter::NORMAL) << typesChanged << " types changed in "
			<< molecule.getName() << typesChanged;
	assert(molecule.getAtom(10).getAtomTypeId()== AtomType::AtomTypeId::P3);
}

/**
 * Simple test- reads in some DHFR ligands and ensures that the same atom types
 *  as are set in the MOL2 input files are discovered
 */
TEST_F(CheckTypesTest, testDhfr) {
	auto molecules = Molecule::readMoleculesFromFile("../resources/dhfr.mol2");
	int moleculeNo = 0;
	for (auto & molecule : molecules) {
		MoleculeSssr moleculeSssr(*molecule);
		moleculeSssr.doSssr();
		molecule->createAtomNeighbourhood();
		int typesChanged = molecule->assignAtomTypes();
		REPORT(Reporter::NORMAL) << typesChanged << " types changed in "
				<< molecule->getName() << typesChanged;
		int nExpected = -1;
		switch (moleculeNo) {
		case 0:
			nExpected = 2;
			break;
		case 1:
			nExpected = 1;
			break;
		case 2:
			nExpected = 0;
			break;
		case 3:
			nExpected = 1;
			break;
		case 4:
			nExpected = 2;
			break;
		case 5:
			nExpected = 1;
			break;
		}
		ASSERT_EQ(typesChanged, nExpected);
		moleculeNo++;
	}
}

}
