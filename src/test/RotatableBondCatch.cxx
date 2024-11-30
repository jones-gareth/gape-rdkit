//
// Created by Gareth Jones on 12/20/2022.
//


#define CATCH_CONFIG_MAIN
#define _USE_MATH_DEFINES

#include "catch2/catch.hpp"
#include "gape/GapeSettings.h"
#include "util/Reporter.h"
#include "mol/MolUtil.h"
#include "gape/RotatableBond.h"
#include <GraphMol/GraphMol.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmilesParse.h>

using namespace Gape;
using namespace Gape;
using namespace RDKit;

SuperpositionMolecule *loadSuperpositionMolecule(const std::string &smilesIn, const GapeSettings &settings) {
    auto mol = SmilesToMol(smilesIn);
    mol->setProp(common_properties::_Name, "Unknown molecule");
    MolOps::addHs(*mol);
    MolOps::sanitizeMol(*mol);
    auto superpositionMolecule = new SuperpositionMolecule(*mol, settings);
    delete mol;
    superpositionMolecule->findFreelyRotatableBonds();
    return superpositionMolecule;
}

std::vector<std::shared_ptr<RotatableBond>>
findRotatableBonds(const std::string &smilesIn, const GapeSettings &settings) {
    auto superpositionMolecule = loadSuperpositionMolecule(smilesIn, settings);
    auto rotatableBonds = superpositionMolecule->getRotatableBonds();
    delete superpositionMolecule;
    return rotatableBonds;
}

TEST_CASE("Find rotatable bonds works as expected", "[RotatableBond]") {
    Reporter::setMinReportingLevel(Reporter::DEBUG);
    const GapeSettings settings;

    SECTION("ETHANE") {
        const auto rotatableBonds = findRotatableBonds("CC", settings);
        CHECK(rotatableBonds.size() == 1);
        CHECK(rotatableBonds[0]->getRotatableBondType() == RotatableBondType::Full);
    }SECTION("PORPANE") {
        const auto rotatableBonds = findRotatableBonds("CCC", settings);
        CHECK(rotatableBonds.size() == 2);
        CHECK(rotatableBonds[0]->getRotatableBondType() == RotatableBondType::Full);
        CHECK(rotatableBonds[1]->getRotatableBondType() == RotatableBondType::Full);
    }SECTION("Amide bond") {
        const auto rotatableBonds = findRotatableBonds("C(=O)NC", settings);
        CHECK(rotatableBonds.size() == 2);
        CHECK(rotatableBonds[0]->getRotatableBondType() == RotatableBondType::Flip);
        CHECK(rotatableBonds[1]->getRotatableBondType() == RotatableBondType::Full);
    }SECTION("PLANAR NITROGEN 1") {
        const auto rotatableBonds = findRotatableBonds("c1(NC)ccccc1", settings);
        CHECK(rotatableBonds.size() == 2);
        CHECK(rotatableBonds[0]->getRotatableBondType() == RotatableBondType::Flip);
        CHECK(rotatableBonds[1]->getRotatableBondType() == RotatableBondType::Full);
    }SECTION("PLANAR NITROGEN 2") {
        const auto rotatableBonds = findRotatableBonds("C=CNC", settings);
        CHECK(rotatableBonds.size() == 2);
        CHECK(rotatableBonds[0]->getRotatableBondType() == RotatableBondType::Flip);
        CHECK(rotatableBonds[1]->getRotatableBondType() == RotatableBondType::Full);
    }SECTION("CARBOXYLIC ACID") {
        const auto rotatableBonds = findRotatableBonds("CC(=O)O", settings);
        CHECK(rotatableBonds.size() == 2);
        CHECK(rotatableBonds[0]->getRotatableBondType() == RotatableBondType::Full);
        CHECK(rotatableBonds[1]->getRotatableBondType() == RotatableBondType::Flip);
    }
}

TEST_CASE("Rotate bonds", "[RotatableBond]") {

    Reporter::setMinReportingLevel(Reporter::DEBUG);
    const GapeSettings settings;

    auto smiles = "O=C1c2cccc3c2[C@H](CCC3)CN1[C@@H]1C[N@H+]2CC[C@@H]1CC2";
    auto superpositionMolecule = loadSuperpositionMolecule(smiles, settings);
    superpositionMolecule->generate3D();
    std::ofstream outFile("rotateBond.sdf");
    outFile << superpositionMolecule->ToMolBlock() << "$$$$" << std::endl;
    auto rotatableBonds = superpositionMolecule->getRotatableBonds();
    auto conformer = superpositionMolecule->getReferenceConformer();
    auto originalConformer = conformer;
    double rotationAngle = M_PI / 2.0;
    auto rotatableBond = rotatableBonds[0];
    rotatableBond->rotateBond(rotationAngle, conformer);
    superpositionMolecule->setConformer(conformer);
    outFile << superpositionMolecule->ToMolBlock() << "$$$$" << std::endl;
    outFile.close();

    const auto tolerance = 1e-5;
    for (const auto &torsion: rotatableBond->getTorsions()) {
        auto expectedAngle = torsion.referenceAngle + rotationAngle;
        const auto actualAngle = torsionAngle(conformer.getAtomPos(torsion.index0),
                                              conformer.getAtomPos(torsion.index1),
                                              conformer.getAtomPos(torsion.index2),
                                              conformer.getAtomPos(torsion.index3));
        if (expectedAngle > M_PI) {
            expectedAngle = -(2.0 * M_PI - expectedAngle);
        }
        if (abs(expectedAngle - actualAngle) > tolerance) {
            std::cerr << "oops" << std::endl;
        }
        CHECK(equals(expectedAngle, actualAngle, tolerance));
    }

    for (int k = 0; k < 2; k++) {
        const auto atomList = k == 0 ? rotatableBond->getAtom1List() : rotatableBond->getAtom2List();
        for (size_t i = 0; i < atomList.size(); i++) {
            for (size_t j = i + 1; j < atomList.size(); j++) {
                const auto idx1 = atomList[i]->getIdx();
                const auto idx2 = atomList[j]->getIdx();
                const auto d1 = (originalConformer.getAtomPos(idx1) - originalConformer.getAtomPos(idx2)).length();
                const auto d2 = (conformer.getAtomPos(idx1) - conformer.getAtomPos(idx2)).length();
                if (abs(d1 - d2) > tolerance) {
                    std::cerr << "oops" << std::endl;
                }
                CHECK(equals(d1, d2, tolerance));
            }
        }
    }
    delete superpositionMolecule;

}
