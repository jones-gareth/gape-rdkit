/*
 * Charge.cpp
 *
 *  Created on: Jan 21, 2016
 *      Author: gjones
 */

#include "Charge.h"
#include "SubstructureSearch.h"
#include "../util/ConfigFile.h"
#include "../util/Reporter.h"

namespace GarethMol {

using namespace std;
using namespace GarethUtil;

/**
 * Reads a line from the charge definition file
 * @param config
 * @return
 */
static unique_ptr<ChargePattern> readConfigFileLine(SmartsParser & parser,
		ConfigFile & config) {
	string pattern = config.readString();
	int formalCharge = config.readInteger();
	double partialCharge = config.readDouble();
	bool multiple = config.readBoolean();

	REPORT(Reporter::DEBUG) << "Adding charge pattern " << pattern
			<< " formal charge " << formalCharge << " partial " << partialCharge
			<< " multiple " << multiple;
	auto chargePattern = make_unique<ChargePattern>(parser, pattern,
			formalCharge, partialCharge, multiple);
	return chargePattern;
}

/**
 * Reads the patterns from the charge file
 *
 * @return
 */
static vector<unique_ptr<ChargePattern>> readChargePatterns() {
	REPORT(Reporter::DEBUG) << "Reading charge patterns ";
	std::vector<unique_ptr<ChargePattern>> vals;
	SmartsParser parser;
	auto func = [& vals, & parser] (ConfigFile & config) {
		vals.push_back(readConfigFileLine(parser, config));
	};
	const string fileName("charge.txt");
	ConfigFile::process(fileName, func);
	return vals;
}

std::vector<unique_ptr<ChargePattern>> & Charge::chargePatterns() {
	static std::vector<unique_ptr<ChargePattern>> patterns =
			readChargePatterns();
	return patterns;
}

void Charge::chargeMolecule(Molecule & molecule) {

	// apply charge patterns
	for (auto & pattern : chargePatterns()) {
		SubstructureSearch ss(*pattern->query, molecule);
		const auto matchFunction =
				[this, &pattern, & molecule] (const std::vector<size_t> & queryIdsToTargetIds) {
					setChargeGroup(*pattern, molecule, queryIdsToTargetIds);
				};
		ss.setCallbackFunction(matchFunction);
		ss.compare();
	}

	// look for other atoms that may need to be charged
	for (auto & atom : molecule.getAtoms()) {

		const auto & nbos = atom->getAtomType().getNeutralBondOrders();
		if (nbos.size() == 0) {
			continue;
		}
		if (atom->getNeighbourhood().isAromatic()) {
			// bond counting is not sufficient for aromatic rings, where we should use
			// patterns or check the ring as a whole
			continue;
		}

		auto currentCharge = atom->getFormalCharge();
		if (currentCharge != 0) {
			REPORT(Reporter::NORMAL) << "Atom " << atom->info()
					<< " currently has formal charge " << currentCharge;
		}
		auto bondOrder = atom->totalBondOrder(molecule);

		auto possibleCharges = mapToNewList<int>(nbos,
				[bondOrder] (int nbo) {return bondOrder-nbo;});
		// not sure now to handle hypervalent atoms- I just look for the smallest
		// differnce bwteen the partial valence and the closest neutral bond order.
		auto charge =
				reduce<int, int>(possibleCharges,
						std::numeric_limits<int>::max(),
						[](int currentCharge, int possibleCharge) {
							return abs(possibleCharge)<abs(currentCharge) ? possibleCharge:currentCharge;
						});
		if (charge != currentCharge) {
			REPORT(Reporter::NORMAL) << "Atom " << atom->info()
					<< " does not satisfy neutral bond order and current charge: setting charge to "
					<< charge;
			atom->setPartialCharge(charge);
			atom->setFormalCharge(charge);
			atom->setFormalChargeSet(charge != 0 ? true : false);
		}
	}
}

void Charge::setChargeGroup(const ChargePattern & pattern, Molecule & molecule,
		const std::vector<size_t> & queryIdsToTargetIds) const {
	// set formal and partial charge on 1st atom
	auto & atom = molecule.getAtom(queryIdsToTargetIds.at(0));
	REPORT(Reporter::NORMAL) << "Atom " << atom.info()
			<< " matches charge pattern " << pattern.pattern
			<< " formal charge " << pattern.formalCharge << " partial charge "
			<< pattern.partialCharge;
	atom.setPartialCharge(pattern.partialCharge);
	if (pattern.formalCharge != 0) {
		atom.setFormalCharge(pattern.formalCharge);
		atom.setFormalChargeSet(true);
	}

	// if the pattern can have only one atom with formal charge set then remove formal charges from other atoms
	if (!pattern.multiple && pattern.formalCharge != 0) {
		for (auto iter = queryIdsToTargetIds.cbegin() + 1;
				iter != queryIdsToTargetIds.cend(); ++iter) {
			auto & atom2 = molecule.getAtom(*iter);
			if (atom2.isFormalChargeSet()) {
				assert(&atom != &atom2);
				REPORT(Reporter::DEBUG) << "Clearing formal charge from "
						<< atom2.info();
				atom2.setFormalCharge(0);
				atom2.setFormalChargeSet(false);
			}
		}
	}
}

} /* namespace GarethMol */
