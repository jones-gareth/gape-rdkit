//
// Created by Gareth Jones on 12/28/2022.
//

#define CATCH_CONFIG_MAIN

#include "catch2/catch.hpp"
#include "gape/GapeApp.h"
#include "util/Reporter.h"
#include "mol/HydrogenBondingType.h"
#include "mol/Solvate.h"
#include <GraphMol/GraphMol.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>

using namespace Gape;
using namespace RDKit;

std::map<Atom *, std::shared_ptr<const HydrogenBondingType>>
findDonorsAndAcceptors(RWMol &mol, const GapeApp &settings, bool solvate = true) {
    mol.setProp(common_properties::_Name, "Unknown molecule");
    MolOps::addHs(mol);
    MolOps::sanitizeMol(mol);
    if (solvate) {
        solvateMolecule(settings.getSolvationRules(), mol);
    }
    auto donorsAndAcceptors = findHydrogenBondDonorsAndAcceptors(settings.getHydrogenBondingTypes(), mol);
    return donorsAndAcceptors;
}

TEST_CASE("Finding hydrogen bonding types works as expected", "[hydrogenBondingType]") {

    Reporter::setMinReportingLevel(Reporter::DEBUG);
    const GapeApp settings;

    SECTION("Terminal Phosphate") {
        auto const mol = "OP(=O)(O)C"_smiles;
        const auto features = findDonorsAndAcceptors(*mol, settings);
        CHECK(features.size() == 3);
        for (const auto &[atom, bondingType]: features) {
            CHECK(atom->getAtomicNum() == 8);
            CHECK(bondingType->hydrogenBondType == HydrogenBondType::Acceptor);
            CHECK(bondingType->name == "Terminal Phosphate");
        }
    }

    SECTION("Phosphinyl") {
        auto const mol = "OP(=O)(C)C"_smiles;
        const auto features = findDonorsAndAcceptors(*mol, settings);
        CHECK(features.size() == 2);
        for (const auto &[atom, bondingType]: features) {
            CHECK(atom->getAtomicNum() == 8);
            CHECK(bondingType->hydrogenBondType == HydrogenBondType::Acceptor);
            CHECK(bondingType->name == "Phosphinyl");
        }
    }

    SECTION("Carboxylate") {
        auto const mol = "OC(=O)C"_smiles;
        const auto features = findDonorsAndAcceptors(*mol, settings);
        CHECK(features.size() == 2);
        for (const auto &[atom, bondingType]: features) {
            CHECK(atom->getAtomicNum() == 8);
            CHECK(bondingType->hydrogenBondType == HydrogenBondType::Acceptor);
            CHECK(bondingType->name == "Carboxylate");
        }
    }

    SECTION("Asymmetric Het6 N") {
        auto const mol = "N1=NC=CC=C1"_smiles;
        const auto features = findDonorsAndAcceptors(*mol, settings);
        CHECK(features.size() == 2);
        for (const auto &[atom, bondingType]: features) {
            CHECK(atom->getAtomicNum() == 7);
            CHECK(bondingType->hydrogenBondType == HydrogenBondType::Acceptor);
            CHECK(bondingType->name == "Asymmetric Het6 N");
        }
    }

    SECTION("Acid O") {
        auto const mol = "O=C(O)C"_smiles;
        const auto features = findDonorsAndAcceptors(*mol, settings, false);
        CHECK(features.size() == 2);
        for (const auto &[atom, bondingType]: features) {
            CHECK(atom->getAtomicNum() == 8);
            CHECK(bondingType->hydrogenBondType == HydrogenBondType::Acceptor);
            if (atom->getTotalNumHs(true) == 1) {
                CHECK(bondingType->name == "Hydroxyl");
            } else {
                CHECK(bondingType->name == "Acid O");
            }
        }
    }

    SECTION("Symmetric Het5 N") {
        auto const mol = "N1=COC=C1"_smiles;
        const auto features = findDonorsAndAcceptors(*mol, settings);
        CHECK(features.size() == 1);
        for (const auto &[atom, bondingType]: features) {
            CHECK(atom->getAtomicNum() == 7);
            CHECK(bondingType->hydrogenBondType == HydrogenBondType::Acceptor);
            CHECK(bondingType->name == "Symmetric Het5 N");
        }
    }

    SECTION("Sulfoxide") {
        auto const mol = "N1=COC=C1"_smiles;
        const auto features = findDonorsAndAcceptors(*mol, settings);
        CHECK(features.size() == 1);
        for (const auto &[atom, bondingType]: features) {
            CHECK(atom->getAtomicNum() == 8);
            CHECK(bondingType->hydrogenBondType == HydrogenBondType::Acceptor);
            CHECK(bondingType->name == "Sulfoxide");
        }
    }

    SECTION("Primary Amine") {
        auto const mol = "N1=COC=C1"_smiles;
        const auto features = findDonorsAndAcceptors(*mol, settings, false);
        CHECK(features.size() == 1);
        for (const auto &[atom, bondingType]: features) {
            CHECK(atom->getAtomicNum() == 7);
            CHECK(bondingType->hydrogenBondType == HydrogenBondType::Acceptor);
            CHECK(bondingType->name == "Primary Amine");
        }
    }

    SECTION("Asymmetric Het5 N") {
        auto const mol = "N1=COC=C1"_smiles;
        const auto features = findDonorsAndAcceptors(*mol, settings);
        CHECK(features.size() == 1);
        for (const auto &[atom, bondingType]: features) {
            CHECK(atom->getAtomicNum() == 7);
            CHECK(bondingType->hydrogenBondType == HydrogenBondType::Acceptor);
            CHECK(bondingType->name == "Asymmetric Het5 N");
        }
    }

    SECTION("Thiocarbonyl") {
        auto const mol = "S=C(C)C"_smiles;
        const auto features = findDonorsAndAcceptors(*mol, settings);
        CHECK(features.size() == 1);
        for (const auto &[atom, bondingType]: features) {
            CHECK(atom->getAtomicNum() == 16);
            CHECK(bondingType->hydrogenBondType == HydrogenBondType::Acceptor);
            CHECK(bondingType->name == "Thiocarbonyl");
        }
    }

    SECTION("Hydroxyl") {
        auto const mol = "OCC"_smiles;
        const auto features = findDonorsAndAcceptors(*mol, settings);
        CHECK(features.size() == 1);
        for (const auto &[atom, bondingType]: features) {
            CHECK(atom->getAtomicNum() == 8);
            CHECK(bondingType->hydrogenBondType == HydrogenBondType::Acceptor);
            CHECK(bondingType->name == "Hydroxyl");
        }
    }
}