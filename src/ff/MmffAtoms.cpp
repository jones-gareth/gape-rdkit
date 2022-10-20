/*
 * MmffAtoms.cpp
 *
 *  Created on: May 15, 2016
 *      Author: gareth
 */

#include "MmffAtoms.h"
#include "../util/ConfigFile.h"
#include "../util/Util.h"

#include <memory>
#include <map>
#include <cassert>

namespace GarethFF {

using namespace GarethUtil;

static unique_ptr<MmffAtomSmarts> readAtomSmartsLine(ConfigFile & config) {
	const string & line = config.getCurrentLine();
	bool ok;
	auto sym = convertString<string>(line, &ok, 0, 3);
	assert(ok);
	trim(sym);
	auto no = convertString<int>(line, &ok, 3, 2);
	assert(ok);
	auto hType = -convertString<int>(line, &ok, 9, 13);
	assert(ok);
	auto mmno = convertString<int>(line, &ok, 6, 2);
	assert(ok);
	auto charge = convertString<int>(line, &ok, 12, 2);
	assert(ok);
	auto valance = convertString<int>(line, &ok, 16, 2);
	assert(ok);
	auto name = convertString<string>(line, &ok, 19, 25);
	assert(ok);
	trim(name);
	auto smi = convertString<string>(line, &ok, 44, 0);
	assert(ok);
	trim(smi);

	auto atomSmarts = make_unique<MmffAtomSmarts>(sym, no, mmno, hType, charge,
			valance, name, smi);
	return atomSmarts;

}

static map<int, unique_ptr<MmffAtomSmarts>> readMmffAtomSmarts() {
	map<int, unique_ptr<MmffAtomSmarts>> atomSmartsMap;
	auto func = [&atomSmartsMap] (ConfigFile & config) {
	    auto smartsPtr = readAtomSmartsLine(config);
		atomSmartsMap.emplace(smartsPtr->getMmff94No(), move(smartsPtr));
	};
}

} /* namespace GarethFF */
