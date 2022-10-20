/*
 * TaffTorsions.cpp
 *
 *  Created on: Mar 15, 2016
 *      Author: gareth
 */

#include "TaffTorsions.h"
#include "Taff.h"

#include "../util/ConfigFile.h"
#include "../util/Reporter.h"

namespace GarethFF {

static unique_ptr<TaffTorsion> readConfigFileLine(ConfigFile & config) {
    auto name1 = config.readString();
    auto name2 = config.readString();
    auto name3 = config.readString();
    auto name4 = config.readString();
    auto bondName = config.readString();

    auto k = config.readDouble();
    auto p = config.readDouble();

    auto atom1TypeId = AtomType::typeFromName(name1).getType();
    auto atom2TypeId = AtomType::typeFromName(name2).getType();
    auto atom3TypeId = AtomType::typeFromName(name3).getType();
    auto atom4TypeId = AtomType::typeFromName(name4).getType();
    auto bondTypeId = BondType::typeFromName(bondName).getType();
    return make_unique<TaffTorsion>(atom1TypeId, atom2TypeId, atom3TypeId,
            atom4TypeId, bondTypeId, k, p);
}

static std::vector<unique_ptr<TaffTorsion>> readTaffTorsions() {

    std::vector<unique_ptr<TaffTorsion>> taffTorsions;
    auto func = [&taffTorsions] (ConfigFile & config) {
        taffTorsions.push_back(readConfigFileLine(config));
    };
    const auto file = "TAFF/TAFF_TORS.txt";
    ConfigFile::process(file, func);
    return taffTorsions;
}

std::vector<unique_ptr<TaffTorsion>> & TaffTorsions::taffTorsionDefinitions() {
    static auto taffTorsionDefinitions = readTaffTorsions();
    return taffTorsionDefinitions;
}

map<int, unique_ptr<MatchedTaffTorsion>> TaffTorsions::parametersForMolecule(
        const Molecule & mol) const {

    map<int, unique_ptr<MatchedTaffTorsion>> matchedTaffTorsionsMap;

    auto nTorsions = mol.getTorsions().size();
    for (auto torsionNo = 0ul; torsionNo < nTorsions; torsionNo++) {
        const auto & torsion = mol.getTorsions().at(torsionNo);

        unique_ptr<MatchedTaffTorsion> matchedTorsion = nullptr;
        auto atom1Type = mol.getAtom(torsion->getAtom1()).getAtomTypeId();
        auto atom2Type = mol.getAtom(torsion->getAtom2()).getAtomTypeId();
        auto atom3Type = mol.getAtom(torsion->getAtom3()).getAtomTypeId();
        auto atom4Type = mol.getAtom(torsion->getAtom4()).getAtomTypeId();
        const auto & bond23 = mol.getBond(torsion->getBond23());

        for (const auto & taffTorsion : taffTorsionDefinitions()) {

            if (!Taff::matchBond(mol, bond23, taffTorsion->getBondTypeId())) {
                continue;
            }

            bool match = false, reverse = false;
            if (Taff::matchAtom(atom1Type, taffTorsion->getAtom1TypeId())
                    && Taff::matchAtom(atom2Type,
                            taffTorsion->getAtom2TypeId())
                    && Taff::matchAtom(atom3Type,
                            taffTorsion->getAtom3TypeId())
                    && Taff::matchAtom(atom4Type,
                            taffTorsion->getAtom4TypeId())) {
                match = true;
            } else if (Taff::matchAtom(atom4Type,
                    taffTorsion->getAtom1TypeId())
                    && Taff::matchAtom(atom3Type,
                            taffTorsion->getAtom2TypeId())
                    && Taff::matchAtom(atom2Type,
                            taffTorsion->getAtom3TypeId())
                    && Taff::matchAtom(atom1Type,
                            taffTorsion->getAtom4TypeId())) {
                match = true;
                reverse = true;
            }

            if (match
                    && (matchedTorsion == nullptr
                            || matchedTorsion->getTaffTorsion()->getWeight()
                                    < taffTorsion->getWeight())) {
                matchedTorsion = make_unique<MatchedTaffTorsion>(reverse,
                        taffTorsion.get(), torsion.get());
            }
        }

        if (matchedTorsion == nullptr) {
            REPORT(Reporter::WARN)
                    << "Unable to find TAFF FF parameters for torsion "
                    << to_string(torsion->getAtom1() + 1) << "["
                    << AtomType::nameFromTypeId(atom1Type) << "] "
                    << to_string(torsion->getAtom2() + 1) << "["
                    << AtomType::nameFromTypeId(atom2Type) << "] "
                    << to_string(torsion->getAtom3() + 1) << "["
                    << AtomType::nameFromTypeId(atom3Type) << "] "
                    << to_string(torsion->getAtom4() + 1) << "["
                    << AtomType::nameFromTypeId(atom4Type) << "]";
        } else {
            matchedTaffTorsionsMap.emplace(torsionNo, move(matchedTorsion));
        }
    }

    return matchedTaffTorsionsMap;
}

const TaffTorsions & TaffTorsions::getInstance() {
    static TaffTorsions taffTorsions;
    return taffTorsions;
}

} /* namespace GarethFF */
