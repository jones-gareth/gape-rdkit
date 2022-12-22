//
// Created by jones on 12/20/2022.
//


#define CATCH_CONFIG_MAIN

#include "catch2/catch.hpp"
#include "gape/GapeApp.h"
#include "util/Reporter.h"
#include "gape/RotatableBond.h"
#include <GraphMol/GraphMol.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmilesParse.h>

using namespace Gape;
using namespace Gape;
using namespace RDKit;

std::vector<std::shared_ptr<const RotatableBond>>
findRotatableBonds(const std::string &smilesIn, const GapeApp &settings) {
    auto mol = SmilesToMol(smilesIn);
    mol->setProp(common_properties::_Name, "Unknown molecule");
    MolOps::addHs(*mol);
    MolOps::sanitizeMol(*mol);
    SuperpositionMolecule superpositionMolecule(*mol, settings);
    delete mol;
    superpositionMolecule.findFreelyRotatableBonds();
    return superpositionMolecule.getRotatableBonds();
}

TEST_CASE("Find rotatable bonds works as expected", "[RotatableBond]") {
    Reporter::setMinReportingLevel(Reporter::DEBUG);
    const GapeApp settings;

    SECTION("ETHANE") {
        const auto rotatableBonds = findRotatableBonds("CC", settings);
        CHECK(rotatableBonds.size() == 1);
        CHECK(rotatableBonds[0]->getRotatableBondType() == RotatableBondType::Full);
    }
    SECTION("PORPANE") {
        const auto rotatableBonds = findRotatableBonds("CCC", settings);
        CHECK(rotatableBonds.size() == 2);
        CHECK(rotatableBonds[0]->getRotatableBondType() == RotatableBondType::Full);
        CHECK(rotatableBonds[1]->getRotatableBondType() == RotatableBondType::Full);
    }
    SECTION("Amide bond") {
        const auto rotatableBonds = findRotatableBonds("C(=O)NC", settings);
        CHECK(rotatableBonds.size() == 2);
        CHECK(rotatableBonds[0]->getRotatableBondType() == RotatableBondType::Flip);
        CHECK(rotatableBonds[1]->getRotatableBondType() == RotatableBondType::Full);
    }
    SECTION("PLANAR NITROGEN 1") {
        const auto rotatableBonds = findRotatableBonds("c1(NC)ccccc1", settings);
        CHECK(rotatableBonds.size() == 2);
        CHECK(rotatableBonds[0]->getRotatableBondType() == RotatableBondType::Flip);
        CHECK(rotatableBonds[1]->getRotatableBondType() == RotatableBondType::Full);
    }
    SECTION("PLANAR NITROGEN 2") {
        const auto rotatableBonds = findRotatableBonds("C=CNC", settings);
        CHECK(rotatableBonds.size() == 2);
        CHECK(rotatableBonds[0]->getRotatableBondType() == RotatableBondType::Flip);
        CHECK(rotatableBonds[1]->getRotatableBondType() == RotatableBondType::Full);
    }
    SECTION("CARBOXYLIC ACID") {
        const auto rotatableBonds = findRotatableBonds("CC(=O)O", settings);
        CHECK(rotatableBonds.size() == 2);
        CHECK(rotatableBonds[0]->getRotatableBondType() == RotatableBondType::Full);
        CHECK(rotatableBonds[1]->getRotatableBondType() == RotatableBondType::Flip);
    }
}