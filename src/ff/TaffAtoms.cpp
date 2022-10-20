/*
 * TaffAtoms.cpp
 *
 *  Created on: Mar 8, 2016
 *      Author: gareth
 */

#include "../util/Reporter.h"
#include "../util/ConfigFile.h"
#include "../mol/Molecule.h"

#include "TaffAtoms.h"
#include "Taff.h"

namespace GarethFF {

using namespace GarethUtil;

static unique_ptr<TaffAtom> readConfigFileLine(ConfigFile & config) {
	auto name = config.readString();
	auto r = config.readDouble();
	auto k = config.readDouble();
	auto atomTypeId = AtomType::typeFromName(name).getType();

	return make_unique<TaffAtom>(r, k, atomTypeId);
}

static std::vector<unique_ptr<TaffAtom>> readTaffAtoms() {

	std::vector<unique_ptr<TaffAtom>> taffAtoms;
	auto func = [&taffAtoms] (ConfigFile & config) {
		taffAtoms.push_back(readConfigFileLine(config));
	};
	const auto file = "TAFF/TAFF_VDW.txt";
	ConfigFile::process(file, func);
	return taffAtoms;
}

std::vector<unique_ptr<TaffAtom>> & TaffAtoms::taffAtomDefinitions()  {
	static auto taffAtomDefinitions = readTaffAtoms();
	return taffAtomDefinitions;
}

map<int, TaffAtom *> TaffAtoms::parametersForMolecule(const Molecule & mol) const{
	map<int, TaffAtom *> taffAtomsMap;
	for (auto i = 0ul; i < mol.nAtoms(); i++) {
		const auto & atom = mol.getAtom(i);
		if (!AtomType::isRealAtom(atom.getAtomTypeId())) {
            continue;
        }
		auto matched = false;
		for (auto & taffAtom : taffAtomDefinitions()) {
			if (Taff::matchAtom(atom.getAtomTypeId(), taffAtom->getAtomTypeId())) {
				taffAtomsMap.emplace(i, taffAtom.get());
				matched = true;
				break;
			}
		}
		if (!matched) {
			REPORT(Reporter::WARN) << "Failed to find definition for atom "
					<< atom.info();
		}
	}
	return taffAtomsMap;
}

const TaffAtoms & TaffAtoms::getInstance() {
	static TaffAtoms taffAtoms;
	return taffAtoms;
}

} /* namespace GarethFF */
