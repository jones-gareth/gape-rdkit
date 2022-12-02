//
// Created by Gareth Jones on 12/1/2022.
//

#define CATCH_CONFIG_MAIN

#include "catch2/catch.hpp"
// #include "/home/packages/Catch2/single_include/catch2/catch.hpp"
#include "gape/GapeApp.h"
#include "util/Reporter.h"
#include "mol/Solvate.h"
#include <GraphMol/GraphMol.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>

using namespace Gape;
using namespace Gape;
using namespace RDKit;

TEST_CASE("Solvation works as expected", "[solvation]") {

    Reporter::setMinReportingLevel(Reporter::DEBUG);
    const GapeApp settings;
    const SolvationRuleList &solvationRules = settings.getSolvationRules();

    SECTION("Primary Aliphatic Amine") {
        std::string smilesIn("CN");
        auto mol = SmilesToMol(smilesIn);
// MolOps::addHs(*mol);
        // MolOps::sanitizeMol(*mol);
        Solvate::solvateMolecule(solvationRules, *mol);
        auto smilesOut = MolToSmiles(*mol);
        REPORT(Reporter::DEBUG) << "Smiles in " << smilesIn << " solvated to " << smilesOut;
        CHECK(smilesOut == "CN[+]");
    }
}