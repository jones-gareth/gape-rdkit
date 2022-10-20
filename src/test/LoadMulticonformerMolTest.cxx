#include <gtest/gtest.h>
#include <stddef.h>
#include <cstdio>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "../mol/Atom.h"
#include "../mol/AtomType.h"
#include "../mol/Bond.h"
#include "../mol/BondType.h"
#include "../mol/MulticonformerMolecule.h"
#include "../mol/Molecule.h"
#include "../util/Reporter.h"
#include "../util/Util.h"

namespace {

using namespace std;
using namespace GarethUtil;
using namespace GarethMol;
using namespace ::testing;

using multiMolPtr = unique_ptr<MulticonformerMolecule>;

/**
 * Test class for reading multiconformer molecules
 */
class LoadMulticonformerMolTest: public Test {

protected:
	LoadMulticonformerMolTest() {
		Reporter::setMinReportingLevel(Reporter::DETAIL);
	}
};

void test1G9v(const vector<multiMolPtr> & mols) {
	ASSERT_EQ(mols.size(), static_cast<size_t>(1));
	auto & mol = mols.at(0);
	mol->initialize();
	ASSERT_EQ(mol->getConformers().size(), static_cast<size_t>(311));

	auto &conformer = mol->getConformer(0);
	CoordVector atom1Coords = conformer.getCoordinates().col(0);
	ASSERT_TRUE(equals(static_cast<double>(atom1Coords[0]), 0.5527, 1));
	ASSERT_TRUE(equals(static_cast<double>(atom1Coords[1]), -4.6763, 1));
	ASSERT_TRUE(equals(static_cast<double>(atom1Coords[2]), 2.8294, 1));
	atom1Coords = conformer.getCoordinates().col(47);
	ASSERT_TRUE(equals(static_cast<double>(atom1Coords[0]), -2.4751, 1));
	ASSERT_TRUE(equals(static_cast<double>(atom1Coords[1]), 2.0890, 1));
	ASSERT_TRUE(equals(static_cast<double>(atom1Coords[2]), 0.7653, 1));

	auto &conformer2 = mol->getConformer(310);
	atom1Coords = conformer2.getCoordinates().col(46);
	ASSERT_TRUE(equals(static_cast<double>(atom1Coords(0)), -0.7781, 1));
	ASSERT_TRUE(equals(static_cast<double>(atom1Coords(1)), -6.7635, 1));
	ASSERT_TRUE(equals(static_cast<double>(atom1Coords[2]), 3.1041, 1));
	atom1Coords = conformer2.getCoordinates().col(47);
	ASSERT_TRUE(equals(static_cast<double>(atom1Coords[0]), -2.9670, 1));
	ASSERT_TRUE(equals(static_cast<double>(atom1Coords[1]), -0.9854, 1));
	ASSERT_TRUE(equals(static_cast<double>(atom1Coords[2]), 4.9868, 1));

}

TEST_F(LoadMulticonformerMolTest, LoadWriteLoad1g9vSdf) {
	const string inFile = "../resources/1g9v_omega.sdf";
	REPORT(Reporter::INFO) << "loading file " << inFile;
	auto mols = MulticonformerMolecule::readMoleculesFromFile(inFile);
	test1G9v(mols);
	const string outFile = "../resources/multiMol.sdf";
	REPORT(Reporter::INFO) << "writing file " << outFile;
	MulticonformerMolecule::writeMoleculesToFile(outFile, mols);
	REPORT(Reporter::INFO) << "loading file " << outFile;
	mols = MulticonformerMolecule::readMoleculesFromFile(outFile);
	test1G9v(mols);
	//remove(outFile.c_str());
}

TEST_F(LoadMulticonformerMolTest, LoadWriteLoad1ywrCcdc) {
	const string inFile = "../resources/1ywr_ccdc.mol2";
	REPORT(Reporter::INFO) << "loading file " << inFile;
	auto mols = MulticonformerMolecule::readMoleculesFromFile(inFile);
	ASSERT_TRUE(mols.size()==1);
}

}
