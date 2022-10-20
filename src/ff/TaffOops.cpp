/*
 * TaffOops.cpp
 *
 *  Created on: Mar 13, 2016
 *      Author: gjones
 */

#include "TaffOops.h"
#include "Taff.h"

#include "../util/Reporter.h"
#include "../util/ConfigFile.h"
#include "../mol/Molecule.h"

namespace GarethFF {

using namespace GarethUtil;

static unique_ptr<TaffOop> readConfigFileLine(ConfigFile & config) {
    auto name = config.readString();
    auto k = config.readDouble();
    auto atomTypeId = AtomType::typeFromName(name).getType();

    return make_unique<TaffOop>(k, atomTypeId);
}

static std::vector<unique_ptr<TaffOop>> readTaffOops() {

    std::vector<unique_ptr<TaffOop>> taffOops;
    auto func = [&taffOops] (ConfigFile & config) {
        taffOops.push_back(readConfigFileLine(config));
    };
    const auto file = "TAFF/TAFF_OOP_BEND.txt";
    ConfigFile::process(file, func);
    return taffOops;
}

std::vector<unique_ptr<TaffOop>> & TaffOops::taffOopDefinitions() {
    static auto taffOopDefinitions = readTaffOops();
    return taffOopDefinitions;
}

map<int, TaffOop *> TaffOops::parametersForMolecule(
        const Molecule & mol) const {

    map<int, TaffOop *> taffAtomsMap;
    for (auto i = 0ul; i < mol.nAtoms(); i++) {
        const auto & atom = mol.getAtom(i);
        if (!AtomType::isRealAtom(atom.getAtomTypeId())) {
            continue;
        }
        for (auto & taffOop : taffOopDefinitions()) {
            if (Taff::matchAtom(atom.getAtomTypeId(), taffOop->getAtomTypeId())) {
                auto nConnections =
                        atom.getNeighbourhood().getAtomNeighbours().size();
                if (nConnections != 3) {
                    stringstream errStream;
                    errStream << "Planar atom "
                            << to_string(atom.getAtomNo() + 1)
                            << " has " << nConnections << " connections!";
                    auto message = errStream.str();
                    if (AtomType::isNitrogenType(atom.getAtomTypeId())) {
                        REPORT(Reporter::DEBUG) << message;
                    } else {
                        REPORT(Reporter::WARN) << message;
                    }
                } else {
                    taffAtomsMap.emplace(i, taffOop.get());
                }
                break;
            }
        }
    }
    return taffAtomsMap;
}

const TaffOops & TaffOops::getInstance() {
    static TaffOops taffOops;
    return taffOops;
}

} /* namespace GarethFF */
