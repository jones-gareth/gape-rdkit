/*
 *SssrTest.cxx
 *
 *	Test file for determining SSSR
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
#include <fstream>

#include "../mol/MoleculeSssr.h"
#include "../mol/Ring.h"
#include "../mol/SmilesParser.h"
#include "../util/Reporter.h"

namespace {

    using namespace std;
    using namespace GarethUtil;
    using namespace GarethMol;
    using namespace ::testing;

    /**
     * Test class for finding SSSR in compounds
     */
    class SssrTest : public Test {
    protected:

        SssrTest() {
            Reporter::setMinReportingLevel(Reporter::NORMAL);
        }
    };

    auto smilesToMolecule(string smiles) {
        SmilesParser parser;
        auto molecule = parser.parseSmiles(smiles);
        MoleculeSssr moleculeSssr(*molecule);
        moleculeSssr.doSssr();
        molecule->createAtomNeighbourhood();
        ofstream out("test.sdf");
        molecule->writeSdfFile(out, "sssr test for " + smiles);
        return molecule;
    }

    TEST_F(SssrTest, testSmarts) {

        string smiles =
                "c1cccc(c1)C23#C4([Co]2([Co]34([CH-]#[O+])([CH-]#[O+])[CH-]#[O+])([CH-]#[O+])([CH-]#[O+])[CH-]#[O+])c5ccccc5";
        auto molecule = smilesToMolecule(smiles);
        auto & rings = molecule->getRings();
        ASSERT_TRUE(rings.size() == 6);
        ASSERT_TRUE(rings.at(0)->getAtoms().size() == 3);
        ASSERT_TRUE(rings.at(1)->getAtoms().size() == 3);
        ASSERT_TRUE(rings.at(2)->getAtoms().size() == 3);
        ASSERT_TRUE(rings.at(3)->getAtoms().size() == 3);
        ASSERT_TRUE(rings.at(4)->getAtoms().size() == 6);
        ASSERT_TRUE(rings.at(5)->getAtoms().size() == 6);

        smiles = "O=S(=O)(CCS(=O)(=O)O)O.[NaH]";
        molecule = smilesToMolecule(smiles);
        ASSERT_TRUE(molecule->getRings().size() == 0);
    }

    //void test2 () {

    TEST_F(SssrTest, testBenzene) {
        const string inFile = "../resources/benzene.mol2";
        Molecule molecule;
        auto ok = molecule.readMoleculeFromFile(inFile);
        ASSERT_TRUE(ok);

        MoleculeSssr moleculeSssr(molecule);
        moleculeSssr.doSssr();
        molecule.createAtomNeighbourhood();

        auto & rings = molecule.getRings();
        ASSERT_TRUE(rings.size() == 1);
        ASSERT_TRUE(rings.at(0)->getAtoms().size() == 6);

    }

    TEST_F(SssrTest, testP0AD64) {
        const string file = "../resources/P0AD64.sdf";
        auto molecules = Molecule::readMoleculesFromFile(file);
        int molNo = 0;
        for (auto & molecule : molecules) {
            MoleculeSssr moleculeSssr(*molecule);
            moleculeSssr.doSssr();
            molecule->createAtomNeighbourhood();
            auto & rings = molecule->getRings();

            switch (molNo) {

                case 0:
                    ASSERT_TRUE(rings.size() == 4);
                    ASSERT_TRUE(rings.at(0)->getAtoms().size() == 5);
                    ASSERT_TRUE(rings.at(1)->getAtoms().size() == 5);
                    ASSERT_TRUE(rings.at(2)->getAtoms().size() == 5);
                    ASSERT_TRUE(rings.at(3)->getAtoms().size() == 7);
                    break;
                case 1:
                    ASSERT_TRUE(rings.size() == 3);
                    ASSERT_TRUE(rings.at(0)->getAtoms().size() == 5);
                    ASSERT_TRUE(rings.at(1)->getAtoms().size() == 6);
                    ASSERT_TRUE(rings.at(2)->getAtoms().size() == 7);
                    break;
                case 2:
                case 3:
                case 4:
                    ASSERT_TRUE(rings.size() == 0);
                    break;
                case 5:
                    ASSERT_TRUE(rings.size() == 1);
                    ASSERT_TRUE(rings.at(0)->getAtoms().size() == 5);
                    break;
                default:
                    throw invalid_argument("switch argument error");
            }
            molNo++;
        }

    }

}
