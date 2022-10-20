/*
 * loadMol.cpp
 *
 *	Test file for reading in molecules
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

#include "../mol/Atom.h"
#include "../mol/AtomType.h"
#include "../mol/Bond.h"
#include "../mol/BondType.h"
#include "../mol/Molecule.h"
#include "../util/Reporter.h"

namespace {

using namespace std;
using namespace GarethUtil;
using namespace GarethMol;
using namespace ::testing;

/**
 * Test class for reading molecules
 */
class LoadMolTest: public Test {

protected:
	LoadMolTest() {
		Reporter::setMinReportingLevel(Reporter::NORMAL);
	}
};

void testHivMols(const vector<Molecule::MoleculePtr> & mols) {
	ASSERT_EQ(mols.size(), static_cast<size_t>(7));

	{
		auto & mol = mols.at(0);
		ASSERT_EQ(mol->nAtoms(), static_cast<size_t>(42));
		ASSERT_EQ(mol->nBonds(), static_cast<size_t>(43));
		ASSERT_EQ(mol->getAtom(0).getAtomType().getType(),
				AtomType::AtomTypeId::C3);
		ASSERT_EQ(mol->getAtom(1).getAtomType().getType(),
				AtomType::AtomTypeId::C2);
		ASSERT_EQ(mol->getAtom(21).getAtomType().getType(),
				AtomType::AtomTypeId::S3);
		ASSERT_EQ(mol->getBond(0).getBondType().getType(),
				BondType::BondTypeId::SINGLE);
		ASSERT_EQ(mol->getBond(1).getBondType().getType(),
				BondType::BondTypeId::AR);
		ASSERT_EQ(mol->getBond(20).getBondType().getType(),
				BondType::BondTypeId::DOUBLE);
		ASSERT_EQ(mol->getBond(12).getBondType().getType(),
				BondType::BondTypeId::AM);

		ASSERT_EQ(mol->getAtom(0).isFormalChargeSet(), false);
		ASSERT_EQ(mol->getAtom(0).isSdfValenceSet(), false);
		ASSERT_EQ(mol->getAtom(0).getSdfValence(), 0);
	}

	{
		auto & mol = mols.at(6);
		ASSERT_EQ(mol->nAtoms(), static_cast<size_t>(33));
		ASSERT_EQ(mol->nBonds(), static_cast<size_t>(35));
		ASSERT_EQ(mol->getAtom(0).getAtomType().getType(),
				AtomType::AtomTypeId::C3);
		ASSERT_EQ(mol->getAtom(1).getAtomType().getType(),
				AtomType::AtomTypeId::F);
		ASSERT_EQ(mol->getAtom(21).getAtomType().getType(),
				AtomType::AtomTypeId::H);
		ASSERT_EQ(mol->getBond(0).getBondType().getType(),
				BondType::BondTypeId::SINGLE);
		ASSERT_EQ(mol->getBond(1).getBondType().getType(),
				BondType::BondTypeId::AR);
		ASSERT_EQ(mol->getBond(20).getBondType().getType(),
				BondType::BondTypeId::DOUBLE);
		ASSERT_EQ(mol->getBond(11).getBondType().getType(),
				BondType::BondTypeId::TRIPLE);
	}
}

TEST_F(LoadMolTest, LoadMol2) {
	const string inFile = "../resources/hiv_rt2.mol2";
	REPORT(Reporter::INFO) << "loading file " << inFile;
	auto mols = Molecule::readMoleculesFromFile(inFile);
	testHivMols(mols);
}

TEST_F(LoadMolTest, LoadSaveLoadMol2) {
	const string inFile = "../resources/hiv_rt2.mol2";
	REPORT(Reporter::INFO) << "loading file " << inFile;
	auto mols = Molecule::readMoleculesFromFile(inFile);
	const string outFile = "../resources/out.mol2";
	REPORT(Reporter::INFO) << "writing file " << outFile;
	Molecule::writeMoleculesToFile(outFile, mols);
	REPORT(Reporter::INFO) << "loading file " << outFile;
	mols = Molecule::readMoleculesFromFile(outFile);
	testHivMols(mols);
	remove(outFile.c_str());
}

void testP0AD64Mols(const vector<Molecule::MoleculePtr>& mols) {
	ASSERT_EQ(mols.size(), static_cast<size_t>(6));

	{
		auto & mol = mols.at(0);
		ASSERT_EQ(mol->nAtoms(), static_cast<size_t>(35));
		ASSERT_EQ(mol->nBonds(), static_cast<size_t>(38));
		ASSERT_EQ(mol->getAtom(0).getAtomType().getType(),
				AtomType::AtomTypeId::S);
		ASSERT_EQ(mol->getAtom(1).getAtomType().getType(),
				AtomType::AtomTypeId::C);
		ASSERT_EQ(mol->getAtom(21).getAtomType().getType(),
				AtomType::AtomTypeId::OANY);
		ASSERT_EQ(mol->getAtom(22).getAtomType().getType(),
				AtomType::AtomTypeId::OANY);
		ASSERT_EQ(mol->getBond(0).getBondType().getType(),
				BondType::BondTypeId::SINGLE);
		ASSERT_EQ(mol->getBond(2).getBondType().getType(),
				BondType::BondTypeId::DOUBLE);
		ASSERT_EQ(mol->getBond(36).getBondType().getType(),
				BondType::BondTypeId::DOUBLE);
		ASSERT_EQ(mol->getBond(37).getBondType().getType(),
				BondType::BondTypeId::SINGLE);
		ASSERT_EQ(mol->getAtom(0).isFormalChargeSet(), false);
		ASSERT_EQ(mol->getAtom(0).isSdfValenceSet(), false);
		ASSERT_EQ(mol->getAtom(0).getSdfValence(), 0);
	}

	{
		auto & mol = mols.at(5);
		ASSERT_EQ(mol->nAtoms(), static_cast<size_t>(33));
		ASSERT_EQ(mol->nBonds(), static_cast<size_t>(33));
		ASSERT_EQ(mol->getAtom(0).getAtomType().getType(),
				AtomType::AtomTypeId::C);
		ASSERT_EQ(mol->getAtom(5).getAtomType().getType(),
				AtomType::AtomTypeId::OANY);
		ASSERT_EQ(mol->getAtom(9).getAtomType().getType(),
				AtomType::AtomTypeId::N);
		ASSERT_EQ(mol->getAtom(22).getAtomType().getType(),
				AtomType::AtomTypeId::H);
		ASSERT_EQ(mol->getBond(0).getBondType().getType(),
				BondType::BondTypeId::SINGLE);
		ASSERT_EQ(mol->getBond(3).getBondType().getType(),
				BondType::BondTypeId::DOUBLE);
		ASSERT_EQ(mol->getBond(23).getBondType().getType(),
				BondType::BondTypeId::DOUBLE);
		ASSERT_EQ(mol->getBond(22).getBondType().getType(),
				BondType::BondTypeId::SINGLE);
	}
}

TEST_F(LoadMolTest, LoadSdf) {
	const string inFile = "../resources/P0AD64.sdf";
	REPORT(Reporter::INFO) << "loading file " << inFile;
	auto mols = Molecule::readMoleculesFromFile(inFile);
	testP0AD64Mols(mols);
}

TEST_F(LoadMolTest, LoadSaveLoadSdf) {
	const string inFile = "../resources/P0AD64.sdf";
	REPORT(Reporter::INFO) << "loading file " << inFile;
	auto mols = Molecule::readMoleculesFromFile(inFile);
	const string outFile = "../resources/out.sdf";
	REPORT(Reporter::INFO) << "writing file " << outFile;
	Molecule::writeMoleculesToFile(outFile, mols);
	REPORT(Reporter::INFO) << "loading file " << outFile;
	mols = Molecule::readMoleculesFromFile(outFile);
	testP0AD64Mols(mols);
	//remove(outFile.c_str());
}

TEST_F(LoadMolTest, Benzene) {
	const string inFile = "../resources/benzene.mol2";
	Molecule molecule;
	auto ok = molecule.readMoleculeFromFile(inFile);
	ASSERT_TRUE(ok);

	molecule.createAtomNeighbourhood();

	const AtomNeighbourhood & neighbourhood =
			molecule.getAtom(2).getNeighbourhood();
	ASSERT_EQ(neighbourhood.countAromaticBonds(), 2);
	ASSERT_EQ(neighbourhood.countSingleAndAmideBonds(), 1);
	ASSERT_EQ(neighbourhood.countDoubleBonds(), 0);
	ASSERT_EQ(neighbourhood.totalBondOrder(), 4);
	ASSERT_EQ(neighbourhood.getAllNeighbours().size(), static_cast<size_t>(3));
	ASSERT_EQ(neighbourhood.getAtomNeighbours().size(), static_cast<size_t>(3));
	ASSERT_EQ(neighbourhood.getBonds().size(), static_cast<size_t>(3));
	ASSERT_EQ(neighbourhood.getAtomBonds().size(), static_cast<size_t>(3));
	ASSERT_EQ(neighbourhood.getHeavyAtomNeighbours().size(),
			static_cast<size_t>(2));
	ASSERT_EQ(neighbourhood.getNoHydrogens(), 1);
}

TEST_F(LoadMolTest, Load1g9vLigand) {
	const string inFile = "../resources/ligand_1g9v.mol";
	Molecule molecule;
	auto ok = molecule.readMoleculeFromFile(inFile);
	ASSERT_TRUE(ok);
	ASSERT_EQ(molecule.nAtoms(), static_cast<size_t>(47));
	ASSERT_EQ(molecule.nBonds(), static_cast<size_t>(48));

	molecule.createAtomNeighbourhood();

	const auto & neighbourhood = molecule.getAtom(1).getNeighbourhood();
	ASSERT_EQ(neighbourhood.countAromaticBonds(), 0);
	ASSERT_EQ(neighbourhood.countSingleAndAmideBonds(), 1);
	ASSERT_EQ(neighbourhood.countDoubleBonds(), 0);
	ASSERT_EQ(neighbourhood.totalBondOrder(), 1);
	ASSERT_EQ(neighbourhood.getAllNeighbours().size(), static_cast<size_t>(1));
	ASSERT_EQ(neighbourhood.getAtomNeighbours().size(), static_cast<size_t>(1));
	ASSERT_EQ(neighbourhood.getBonds().size(), static_cast<size_t>(1));
	ASSERT_EQ(neighbourhood.getAtomBonds().size(), static_cast<size_t>(1));
	ASSERT_EQ(neighbourhood.getHeavyAtomNeighbours().size(),
			static_cast<size_t>(1));
	ASSERT_EQ(neighbourhood.getNoHydrogens(), 0);

	const auto & neighbourhood2 = molecule.getAtom(2).getNeighbourhood();
	ASSERT_EQ(neighbourhood2.countAromaticBonds(), 0);
	ASSERT_EQ(neighbourhood2.countSingleAndAmideBonds(), 0);
	ASSERT_EQ(neighbourhood2.countDoubleBonds(), 1);
	ASSERT_EQ(neighbourhood2.totalBondOrder(), 2);
	ASSERT_EQ(neighbourhood2.getAllNeighbours().size(), static_cast<size_t>(1));
	ASSERT_EQ(neighbourhood2.getAtomNeighbours().size(),
			static_cast<size_t>(1));
	ASSERT_EQ(neighbourhood2.getBonds().size(), static_cast<size_t>(1));
	ASSERT_EQ(neighbourhood2.getAtomBonds().size(), static_cast<size_t>(1));
	ASSERT_EQ(neighbourhood2.getHeavyAtomNeighbours().size(),
			static_cast<size_t>(1));
	ASSERT_EQ(neighbourhood2.getNoHydrogens(), 0);

	const auto & neighbourhood3 = molecule.getAtom(28).getNeighbourhood();
	ASSERT_EQ(neighbourhood3.countAromaticBonds(), 0);
	ASSERT_EQ(neighbourhood3.countSingleAndAmideBonds(), 3);
	ASSERT_EQ(neighbourhood3.countDoubleBonds(), 0);
	ASSERT_EQ(neighbourhood3.totalBondOrder(), 3);
	ASSERT_EQ(neighbourhood3.getAllNeighbours().size(), static_cast<size_t>(3));
	ASSERT_EQ(neighbourhood3.getAtomNeighbours().size(),
			static_cast<size_t>(3));
	ASSERT_EQ(neighbourhood3.getBonds().size(), static_cast<size_t>(3));
	ASSERT_EQ(neighbourhood3.getAtomBonds().size(), static_cast<size_t>(3));
	ASSERT_EQ(neighbourhood3.getHeavyAtomNeighbours().size(),
			static_cast<size_t>(2));
	ASSERT_EQ(neighbourhood3.getNoHydrogens(), 1);

	const auto & neighbourhood4 = molecule.getAtom(38).getNeighbourhood();
	ASSERT_EQ(neighbourhood4.countAromaticBonds(), 0);
	ASSERT_EQ(neighbourhood4.countSingleAndAmideBonds(), 2);
	ASSERT_EQ(neighbourhood4.countDoubleBonds(), 1);
	ASSERT_EQ(neighbourhood4.totalBondOrder(), 4);
	ASSERT_EQ(neighbourhood4.getAllNeighbours().size(), static_cast<size_t>(3));
	ASSERT_EQ(neighbourhood4.getAtomNeighbours().size(),
			static_cast<size_t>(3));
	ASSERT_EQ(neighbourhood4.getBonds().size(), static_cast<size_t>(3));
	ASSERT_EQ(neighbourhood4.getAtomBonds().size(), static_cast<size_t>(3));
	ASSERT_EQ(neighbourhood4.getHeavyAtomNeighbours().size(),
			static_cast<size_t>(3));
	ASSERT_EQ(neighbourhood4.getNoHydrogens(), 0);

}

}
