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

std::map<Atom*, std::shared_ptr<const HydrogenBondingType>> findDonorsAndAcceptors(RWMol mol, const GapeApp &settings) {
    mol.setProp(common_properties::_Name, "Unknown molecule");
    MolOps::addHs(mol);
    MolOps::sanitizeMol(mol);
    solvateMolecule(settings.getSolvationRules(), mol);
    auto donorsAndAcceptors = findHydrogenBondDonorsAndAcceptors(settings.getHydrogenBondingTypes(), mol);
    return donorsAndAcceptors;
}

TEST_CASE("Finding hydrogen bonding types works as expected", "[hydrogenBondingType]") {

    Reporter::setMinReportingLevel(Reporter::DEBUG);
    const GapeApp settings;

    SECTION("Terminal Phosphate") {
        auto const mol = "OP(=O)(O)C"_smiles;
        const auto features = findDonorsAndAcceptors(*mol, settings);
        for (const auto &[atom, bondingType]: features) {
            CHECK(atom->getAtomicNum() == 8);
            CHECK(bondingType->hydrogenBondType == HydrogenBondType::Acceptor);
            CHECK(bondingType->name == "Terminal Phosphate");
        }
    }


}