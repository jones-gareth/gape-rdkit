/*
 * TaffAngle.cpp
 *
 *  Created on: Mar 12, 2016
 *      Author: gjones
 */

#include "../util/Reporter.h"
#include "../util/ConfigFile.h"
#include "../mol/Molecule.h"
#include "TaffAngles.h"
#include "Taff.h"

namespace GarethFF {

using namespace GarethUtil;
using AtomTypeId = AtomType::AtomTypeId;

static unique_ptr<TaffAngle> readConfigFileLine(ConfigFile & config) {
    auto name1 = config.readString();
    auto name2 = config.readString();
    auto name3 = config.readString();
    auto angle = config.readDouble();
    auto k = config.readDouble();
    auto atom1TypeId = AtomType::typeFromName(name1).getType();
    auto atom2TypeId = AtomType::typeFromName(name2).getType();
    auto atom3TypeId = AtomType::typeFromName(name3).getType();

    return make_unique<TaffAngle>(atom1TypeId, atom2TypeId, atom3TypeId, angle,
            k);
}

static std::vector<unique_ptr<TaffAngle>> readTaffAngles() {
    std::vector<unique_ptr<TaffAngle>> taffAngles;
    auto func = [&taffAngles] (ConfigFile & config) {
        taffAngles.push_back(readConfigFileLine(config));
    };
    const auto file = "TAFF/TAFF_ANGLE_BEND.txt";
    ConfigFile::process(file, func);
    return taffAngles;
}

std::vector<unique_ptr<TaffAngle>> & TaffAngles::taffAngleDefinitions() {
    static auto taffAngleDefinitions = readTaffAngles();
    return taffAngleDefinitions;
}

map<int, TaffAngle *> TaffAngles::parametersForMolecule(
        const Molecule & mol) const {
    auto & taffAngles = taffAngleDefinitions();
    map<int, TaffAngle *> taffAnglesMap;

    const auto & angles = mol.getAngles();
    for (auto angleNo = 0ul; angleNo < angles.size(); angleNo++) {
        const auto & angle = angles.at(angleNo);

        TaffAngle * matchedAngle = nullptr;
        auto atomType1 = mol.getAtom(angle->getAtom1()).getAtomTypeId();
        auto atomType2 = mol.getAtom(angle->getAtom2()).getAtomTypeId();
        auto atomType3 = mol.getAtom(angle->getAtom3()).getAtomTypeId();

        if (!AtomType::isRealAtom(atomType1)) {
            continue;
        }
        if (!AtomType::isRealAtom(atomType2)) {
            continue;
        }
        if (!AtomType::isRealAtom(atomType3)) {
            continue;
        }

        for (const auto & taffAngle : taffAngles) {

            if (!Taff::matchAtom(atomType2, taffAngle->getAtom2TypeId())) {
                continue;
            }

            TaffAngle * testAngle = nullptr;
            if (Taff::matchAtom(atomType1, taffAngle->getAtom1TypeId())
                    && Taff::matchAtom(atomType3, taffAngle->getAtom3TypeId())) {
                testAngle = taffAngle.get();
            }

            if (Taff::matchAtom(atomType3, taffAngle->getAtom1TypeId())
                    && Taff::matchAtom(atomType1, taffAngle->getAtom3TypeId())) {
                testAngle = taffAngle.get();
            }

            if (testAngle != nullptr
                    && (matchedAngle == nullptr
                            || matchedAngle->getWeight()
                                    < testAngle->getWeight())) {
                matchedAngle = testAngle;
            }
        }

        if (matchedAngle == nullptr) {
            matchedAngle =
                    defaultParameters().at(
                            mol.getAtom(angle->getAtom2()).getAtomType().getGeometry()).get();
            REPORT(Reporter::WARN) << "Failed to match angle between atoms "
                    << to_string(angle->getAtom1() + 1) << ", "
                    << to_string(angle->getAtom2() + 1) << " and "
                    << to_string(angle->getAtom3() + 1) << " using angle "
                    << matchedAngle->getAngle() << " k "
                    << matchedAngle->getK();
        }
        assert(matchedAngle != nullptr);
        taffAnglesMap.emplace(angleNo, matchedAngle);
    }

    return taffAnglesMap;
}

map<AtomType::Geometry, unique_ptr<TaffAngle>> createDefaultGeometries() {
    map<AtomType::Geometry, unique_ptr<TaffAngle>> geometryMap;

    geometryMap.emplace(AtomType::Geometry::LIN,
            make_unique<TaffAngle>(AtomTypeId::ANY, AtomTypeId::ANY,
                    AtomTypeId::ANY, 180, 0.04));
    geometryMap.emplace(AtomType::Geometry::TRI,
            make_unique<TaffAngle>(AtomTypeId::ANY, AtomTypeId::ANY,
                    AtomTypeId::ANY, 120, 0.03));
    geometryMap.emplace(AtomType::Geometry::TET,
            make_unique<TaffAngle>(AtomTypeId::ANY, AtomTypeId::ANY,
                    AtomTypeId::ANY, 109.5, 0.02));

    return geometryMap;
}

map<AtomType::Geometry, unique_ptr<TaffAngle>> & TaffAngles::defaultParameters() {
    static auto geometryMap = createDefaultGeometries();
    return geometryMap;
}

TaffAngles & TaffAngles::getInstance() {
    static TaffAngles taffAngles;
    return taffAngles;
}

} /* namespace GarethFF */
