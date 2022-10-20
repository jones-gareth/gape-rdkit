/*
 * ChargeMolTest.cpp
 *
 *  Created on: Sep 10, 2015
 *      Author: gjones
 */

#include <gtest/gtest.h>
#include "../util/Reporter.h"
#include "../mol/Molecule.h"
#include "../mol/Charge.h"

namespace {

using namespace std;
using namespace GarethUtil;
using namespace GarethMol;
using namespace ::testing;

/**
 * Tests charging of molecules
 */
class ChargeMolTest: public Test {

protected:

	static void SetUpTestCase() {
		Reporter::setMinReportingLevel(Reporter::DEBUG);
		//Reporter::setFileStream(out);
	}

	static void TearDownTestCase() {
		//out.close();
	}

private:
	static ofstream out;
};

ofstream ChargeMolTest::out { "ChargeMolTest.log" };

unique_ptr<Molecule> processFile(string fileName) {
	REPORT(Reporter::NORMAL) << "Processing file " << fileName;
	auto molecule = Molecule::readAndInitializeMolecule(fileName);
	Charge charge;
	charge.chargeMolecule(*molecule);
	return molecule;
}

TEST_F(ChargeMolTest, Test1ia1) {
	auto molecule = processFile("../resources/ligand_1ia1.mol");

	for (const auto & atom : molecule->getAtoms()) {
		auto formalCharge = atom->getFormalCharge();
		auto partialCharge = atom->getPartialCharge();
		REPORT(Reporter::DEBUG) << "Atom " << atom->info() << " formal charge "
				<< formalCharge << " partial charge " << partialCharge;
		if (atom->getAtomNo() == 7) {
			ASSERT_EQ(formalCharge, 0);
			ASSERT_EQ(partialCharge, 0.5);
		} else if (atom->getAtomNo() == 5) {
			ASSERT_EQ(formalCharge, 1);
			ASSERT_EQ(partialCharge, 0.5);
		} else {
			ASSERT_EQ(formalCharge, 0);
			ASSERT_EQ(partialCharge, 0);
		}
	}
}

TEST_F(ChargeMolTest, Test1jd0) {
	auto molecule = processFile("../resources/ligand_1jd0.mol");

	for (const auto & atom : molecule->getAtoms()) {
		auto formalCharge = atom->getFormalCharge();
		auto partialCharge = atom->getPartialCharge();
		REPORT(Reporter::DEBUG) << "Atom " << atom->info() << " formal charge "
				<< formalCharge << " partial charge " << partialCharge;
		if (atom->getAtomNo() == 7) {
			ASSERT_EQ(formalCharge, -1);
			ASSERT_EQ(partialCharge, -1);
		} else {
			ASSERT_EQ(formalCharge, 0);
			ASSERT_EQ(partialCharge, 0);
		}
	}
}

TEST_F(ChargeMolTest, Test1uou) {
	auto molecule = processFile("../resources/ligand_1uou.mol");

	for (const auto & atom : molecule->getAtoms()) {
		auto formalCharge = atom->getFormalCharge();
		auto partialCharge = atom->getPartialCharge();
		REPORT(Reporter::DEBUG) << "Atom " << atom->info() << " formal charge "
				<< formalCharge << " partial charge " << partialCharge;
		if (atom->getAtomNo() == 17) {
			ASSERT_EQ(formalCharge, 1);
			ASSERT_EQ(partialCharge, 1);
		} else if (atom->getAtomNo() == 20) {
			ASSERT_EQ(formalCharge, -1);
			ASSERT_EQ(partialCharge, -1);
		} else {
			ASSERT_EQ(formalCharge, 0);
			ASSERT_EQ(partialCharge, 0);
		}
	}
}

TEST_F(ChargeMolTest, Test1r58) {
	auto molecule = processFile("../resources/ligand_1r58.mol");

	for (const auto & atom : molecule->getAtoms()) {
		auto formalCharge = atom->getFormalCharge();
		auto partialCharge = atom->getPartialCharge();
		REPORT(Reporter::DEBUG) << "Atom " << atom->info() << " formal charge "
				<< formalCharge << " partial charge " << partialCharge;
		ASSERT_EQ(formalCharge, 0);
		ASSERT_EQ(partialCharge, 0);
	}

}

TEST_F(ChargeMolTest, Test1l7f) {
	auto molecule = processFile("../resources/ligand_1l7f.mol");

	for (const auto & atom : molecule->getAtoms()) {
		auto formalCharge = atom->getFormalCharge();
		auto partialCharge = atom->getPartialCharge();
		REPORT(Reporter::DEBUG) << "Atom " << atom->info() << " formal charge "
				<< formalCharge << " partial charge " << partialCharge;
		if (atom->getAtomNo() == 3) {
			ASSERT_EQ(formalCharge, 1);
			ASSERT_EQ(partialCharge, 0.33333);
		} else if (atom->getAtomNo() == 1) {
			ASSERT_EQ(formalCharge, 0);
			ASSERT_EQ(partialCharge, 0.33333);
		} else if (atom->getAtomNo() == 4) {
			ASSERT_EQ(formalCharge, 0);
			ASSERT_EQ(partialCharge, 0.33333);
		} else if (atom->getAtomNo() == 20) {
			ASSERT_EQ(formalCharge, 0);
			ASSERT_EQ(partialCharge, -0.5);
		} else if (atom->getAtomNo() == 21) {
			ASSERT_EQ(formalCharge, -1);
			ASSERT_EQ(partialCharge, -0.5);
		} else {
			ASSERT_EQ(formalCharge, 0);
			ASSERT_EQ(partialCharge, 0);
		}
	}
}

}
