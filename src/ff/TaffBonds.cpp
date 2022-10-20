/*
 * TaffBonds.cpp
 *
 *  Created on: Mar 12, 2016
 *      Author: gjones
 */

#include "../util/Reporter.h"
#include "../util/ConfigFile.h"
#include "TaffBonds.h"
#include "Taff.h"

namespace GarethFF {

using namespace GarethUtil;
using AtomTypeId = AtomType::AtomTypeId;

static unique_ptr<TaffBond> readConfigFileLine(ConfigFile & config) {
    auto name1 = config.readString();
    auto name2 = config.readString();
    auto bondName = config.readString();
    auto length = config.readDouble();
    auto k = config.readDouble();
    auto atom1TypeId = AtomType::typeFromName(name1).getType();
    auto atom2TypeId = AtomType::typeFromName(name2).getType();
    auto bondType = BondType::typeFromName(bondName).getType();

    return make_unique<TaffBond>(atom1TypeId, atom2TypeId, bondType, length, k);
}

static std::vector<unique_ptr<TaffBond>> readTaffBonds() {
    std::vector<unique_ptr<TaffBond>> taffBonds;
    auto func = [&taffBonds] (ConfigFile & config) {
        taffBonds.push_back(readConfigFileLine(config));
    };
    const auto file = "TAFF/TAFF_BOND_STRETCH.txt";
    ConfigFile::process(file, func);
    return taffBonds;
}

std::vector<unique_ptr<TaffBond>> & TaffBonds::taffBondDefinitions() {
    static auto taffBondDefinitions = readTaffBonds();
    return taffBondDefinitions;
}


map<int, TaffBond *> TaffBonds::parametersForMolecule(
        const Molecule & mol) const {

    auto & taffBonds = taffBondDefinitions();
    map<int, TaffBond *> taffBondsMap;

    for (auto bondNo = 0ul; bondNo < mol.nBonds(); bondNo++) {
        const auto & bond = mol.getBond(bondNo);
        const auto atom1Type = bond.getAtom1().getAtomTypeId();
        const auto atom2Type = bond.getAtom2().getAtomTypeId();

        if (!AtomType::isRealAtom(atom1Type)) {
            continue;
        }
        if (!AtomType::isRealAtom(atom2Type)) {
            continue;
        }

        TaffBond * match = nullptr;
        for (const auto & taffBond : taffBonds) {

            if (!Taff::matchBond(mol, bond, taffBond->getBondTypeId())) {
                continue;
            }

            TaffBond * testBond = nullptr;
            if (Taff::matchAtom(atom1Type, taffBond->getAtom1TypeId())
                    && Taff::matchAtom(atom2Type, taffBond->getAtom2TypeId())) {
                testBond = taffBond.get();
            }
            if (Taff::matchAtom(atom2Type, taffBond->getAtom1TypeId())
                    && Taff::matchAtom(atom1Type, taffBond->getAtom2TypeId())) {
                testBond = taffBond.get();
            }

            if (testBond != nullptr
                    && (match == nullptr
                            || match->getWeight() < testBond->getWeight())) {
                match = testBond;
            }
        }

        if (match == nullptr) {
            match = defaultParameters().get();
            REPORT(Reporter::WARN)
                    << "Unable to find taff bond parameters to match bond between atoms "
                    << to_string(bond.getAtom1().getAtomNo() + 1) << " and "
                    << to_string(bond.getAtom2().getAtomNo() + 1)
                    << " using length " << match->getLength() << " k "
                    << match->getK();
        }
        assert(match != nullptr);
        taffBondsMap.emplace(bondNo, match);
    }

    return taffBondsMap;
}

unique_ptr<TaffBond> & TaffBonds::defaultParameters() {
    static auto defaultParameters = make_unique<TaffBond>(AtomTypeId::WILD,
            AtomTypeId::WILD, BondType::BondTypeId::ANY, 1.5, 600);
    return defaultParameters;
}

const TaffBonds & TaffBonds::getInstance() {
    static TaffBonds taffBonds;
    return taffBonds;
}
} /* namespace GarethFF */
