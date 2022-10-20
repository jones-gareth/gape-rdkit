/*
 * UllmanIsomorphismTest.cxx
 *
 *  Created on: Sep 8, 2015
 *      Author: gjones
 */
#include <gtest/gtest.h>
#include "../mol/Molecule.h"
#include "../mol/MolComparer.h"

namespace {

using namespace std;
using namespace GarethUtil;
using namespace GarethMol;
using namespace ::testing;

/**
 * Test class for subgraph isomorphisms
 */
class UllmanIsomorphismTest: public Test {

protected:
	UllmanIsomorphismTest() {
		Reporter::setMinReportingLevel(Reporter::DETAIL);
	}
};

/**
 * Test of some hiv inhibitors.  Expected values obtained from Java code
 */
TEST_F(UllmanIsomorphismTest, CompareHiv) {
	const string inFile = "../resources/hiv_rt2.mol2";
	const map<string, size_t> results { { "pdb1bqm", 2 }, { "pdb1klm", 8 }, {
			"pdb1rt3", 2 }, { "pdb1vru", 2 }, { "pdb1tvr", 2 },
			{ "pdb1dtt", 1 }, { "pdb1fk9", 12 } };
	auto molecules = Molecule::readMoleculesFromFile(inFile);
	for (auto & mol : molecules) {
		REPORT(Reporter::DETAIL) << "Processing molecule " << mol->getName();
		auto comparer = MolComparer::createIsomorphismComparer(*mol);
		auto nIsomorphisms = comparer->compare();
		ASSERT_EQ(results.at(mol->getName()), nIsomorphisms);
	}
}

/**
 * Test of benzene
 */
TEST_F(UllmanIsomorphismTest, CompareBenzene) {
	const string inFile = "../resources/benzene.mol2";

	auto molecules = Molecule::readMoleculesFromFile(inFile);
	assert(molecules.size() == 1);
	REPORT(Reporter::DETAIL) << "Processing Benzene ";
	auto comparer = MolComparer::createIsomorphismComparer(*molecules.at(0));
	auto nIsomorphisms = comparer->compare();
	ASSERT_EQ(nIsomorphisms, static_cast<size_t>(12));
	ASSERT_EQ(comparer->getQueryAtoms().size(), static_cast<size_t>(6));
	ASSERT_EQ(comparer->getTargetAtoms().size(), static_cast<size_t>(6));
}

}
