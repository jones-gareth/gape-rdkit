//
// Created by Gareth Jones on 12/1/2022.
//

#define CATCH_CONFIG_MAIN

#include "catch2/catch.hpp"
#include "../gape/GapeApp.h"
#include "../util/Reporter.h"
#include "mol/Solvate.h"
#include <GraphMol/GraphMol.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>

using namespace Gape;
using namespace RDKit;

void
solvateSmiles(const std::string &smilesIn, const std::string expectedSmiles, const SolvationRuleList &solvationRules) {
    auto mol = SmilesToMol(smilesIn);
    mol->setProp(common_properties::_Name, "Unknown molecule");
    MolOps::addHs(*mol);
    MolOps::sanitizeMol(*mol);
    solvateMolecule(solvationRules, *mol);
    MolOps::removeAllHs(*mol);
    auto smilesOut = MolToSmiles(*mol);

    auto expectedMol = SmilesToMol(expectedSmiles);
    auto expected = MolToSmiles(*expectedMol);

    REPORT(Reporter::DEBUG) << "Smiles " << smilesIn << " solvated to " << smilesOut << " Expected " << expected;
    CHECK(smilesOut == expected);
    delete mol;
    delete expectedMol;
}

TEST_CASE("Solvation works as expected", "[solvation]") {

    Reporter::setMinReportingLevel(Reporter::DEBUG);
    const GapeApp settings;
    const SolvationRuleList &solvationRules = settings.getSolvationRules();

    SECTION("Primary Aliphatic Amine") {
        solvateSmiles("CN", "[NH3+]C", solvationRules);
    }

    SECTION("Secondary Aliphatic Amine") {
        solvateSmiles("CNC", "C[NH2+]C", solvationRules);
    }

    SECTION("Tertiary Aliphatic Amine") {
        solvateSmiles("CN(C)C", "C[NH+](C)C", solvationRules);
    }

    SECTION("Guanidinium") {
        solvateSmiles("N=C(N)N", "[NH2+]=C(N)(N)", solvationRules);
    }

    SECTION("Amidine 1") {
        solvateSmiles("[NH]=C(N)C", "[NH2+]=C(N)C", solvationRules);
    }

    SECTION("Amidine 2") {
        solvateSmiles("CN=C(N)C", "C[NH+]=C(N)C", solvationRules);
    }

    SECTION("Ortho Aminopryidine") {
        solvateSmiles("n1c(N)cncc1", "[nH+]1c(N)cncc1", solvationRules);
    }

    SECTION("Ortho Aminopryidine 2") {
        solvateSmiles("n1c(N)cccc1", "[nH+]1c(N)cccc1", solvationRules);
    }

    SECTION("Para Aminopryidine") {
        solvateSmiles("n1cnc(-N)cc1", "[nH+]1cnc(N)cc1", solvationRules);
    }

    SECTION("Para Aminopryidine 2") {
        solvateSmiles("n1ccc(-N)cc1", "[nH+]1ccc(N)cc1", solvationRules);
    }

    SECTION("Carboxyl") {
        solvateSmiles("OC=O", "[O-]C=O", solvationRules);
    }

    SECTION("Nitro") {
        solvateSmiles("ON=O", "[O-]N=O", solvationRules);
    }

    SECTION("Tetrazole Tautomer 1") {
        solvateSmiles("[nH]1nncn1", "[n-]1nncn1", solvationRules);
    }

    SECTION("Tetrazole Tautomer 2") {
        solvateSmiles("[nH]1nnnc1", "[n-]1nnnc1", solvationRules);
    }

    SECTION("Protonated Purine") {
        solvateSmiles("[nH]1cnc2c1ncnc2", "[n-]1cnc2c1ncnc2", solvationRules);
    }

    SECTION("Acylsulponamide") {
        solvateSmiles("N(C=O)S(=O)=O", "[N-](C=O)S(=O)=O", solvationRules);
    }

    SECTION("Wyeth heterocycle") {
        solvateSmiles("N1C(=O)ONC(=O)1", "[N-]1C(=O)ONC(=O)1", solvationRules);
    }

    SECTION("Phosphate") {
        solvateSmiles("OP=O", "[O-]P=O", solvationRules);
    }

    SECTION("Sulphonic acid") {
        solvateSmiles("OS(=O)=O", "[O-]S(=O)=O", solvationRules);
    }
}