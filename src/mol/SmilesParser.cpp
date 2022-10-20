/*
 * SmilesParser.cpp
 *
 *  Created on: Sep 22, 2015
 *      Author: gjones
 */

#include "SmilesParser.h"

#include <boost/format.hpp>
#include <cassert>
#include <cctype>
#include <initializer_list>
#include <iostream>
#include <stdexcept>
#include <utility>
#include <vector>

#include "../util/Reporter.h"
#include "Atom.h"
#include "AtomType.h"
#include "Bond.h"
#include "BondType.h"

namespace GarethMol {

using namespace std;
using namespace GarethUtil;

SmilesParser::SmilesParser() {
}

SmilesParser::~SmilesParser() {
}

unique_ptr<Molecule> SmilesParser::parseSmiles(const string & smiles_) {

	auto pos = smiles_.find(' ');
	string name = "Smiles";
	smiles = pos == string::npos ? smiles_ : smiles_.substr(0, pos);
	if (pos != string::npos) {
		name = smiles_.substr(pos + 1, string::npos);
	}

	position = 0;
	atomNo = 0;
	bondNo = 0;
	atoms.clear();
	bonds.clear();
	ringOpenings.clear();
	openBranches = {};
	atoms.reserve(smiles.length() * 2);
	bonds.reserve(smiles.length() * 2);
	currentAtomAromatic = false;
	aromaticAtomNumbers = {};
	currentAtom = nullptr;
	currentBond = BondType::BondTypeId::UNK;
	firstAtomRead = false;

	// parse atom entries
	while (position < smiles.length()) {
		parseAtomEntry();
	}
	// add implicit hydrogens
	addHydrogens();

	auto molecule = make_unique<Molecule>(name + " " + smiles, atoms, bonds);
	return molecule;
}

void SmilesParser::addHydrogens() {
	// iterate over array of  underlying pointer type as if atom array is reallocated
	// as hydrogens are added then unique_pointer will be invalid.
    auto map = [] (const unique_ptr<Atom> & atom) -> Atom * {
        return atom.get();
    };
    auto atomPtrs = mapToNewList<unique_ptr<Atom>, Atom *>(atoms, map);
    for (auto atom: atomPtrs) {
		auto bondOrder = atom->totalBondOrder(bonds);
		auto partialValance = bondOrder - atom->getFormalCharge();
		const auto & nbos = atom->getAtomType().getNeutralBondOrders();
		if (nbos.size() == 0) {
			continue;
		}
		int nHydrogens = -1;
		int nbo = -1;
		for (auto iter = nbos.cbegin(); iter != nbos.cend(); ++iter) {
			nbo = *iter;
			if (nbo >= partialValance) {
				nHydrogens = nbo - partialValance;
				break;
			}
		}

		if (nHydrogens < 0) {
			auto message =
					parseMessage(
							(boost::format(
									"Valence is too high for atom %s: bond order %d charge %d partial valence %d NBOS [%s] NBO %d nHydrogens %d")
									% atom->info() % bondOrder
									% atom->getFormalCharge() % partialValance
									% collectionToString(nbos, ", ") % nbo
									% nHydrogens).str());
			REPORT(Reporter::WARN) << message;
			nHydrogens = 0;
		}
		if (nHydrogens > 0) {
			addHydrogens(atom, nHydrogens, true);
			REPORT(Reporter::DEBUG)
					<< boost::format("added %d implicit hydrogens to atom %s")
							% nHydrogens % atom->info();
		}
	}
}

void SmilesParser::parseAtomEntry() {
	if (firstAtomRead) {
		// deal with disconnected structures
		if (getCurrentCharacter() == '.') {
			position++;
			currentAtom = nullptr;
			if (currentBond != BondType::BondTypeId::UNK) {
				parseError("Unused bond at structure break!");
			}
		} else {
			// parse bond and open branch prior to next atom
			parseOpenBranch();
			if (currentBond == BondType::BondTypeId::UNK) {
				parseBond();
			}
		}
	}
	parseAtom();
	do {
		parseBond();
	} while (parseRingClosure());
	parseCloseBranch();
}

void SmilesParser::parseCloseBranch() {
	// handle closed branch
	while (getCurrentCharacter() == ')') {
		position++;
		if (openBranches.size() == 0) {
			parseError("No matching open branch for close");
		}
		REPORT(Reporter::DEBUG) << parseMessage("parsed close branch");
		// pop current atom from stack of open branches
		currentAtom = openBranches.back();
		openBranches.pop_back();
	}
}

void SmilesParser::parseOpenBranch() {
	// handle open branch
	if (getCurrentCharacter() == '(') {
		position++;
		// add current atom to stack of open branches
		openBranches.push_back(currentAtom);
		REPORT(Reporter::DEBUG) << parseMessage("parsed open branch");
	}
}

bool SmilesParser::parseRingClosure() {
	// handle ring opeening and closure

	auto currentChar = getCurrentCharacter();
	if (isdigit(currentChar) || currentChar == '%') {
		// find ring number
		int ringClosureNo = -1;
		if (isdigit(currentChar)) {
			ringClosureNo = currentChar - '0';
			position++;
		} else if (currentChar == '%') {
			position++;
			auto char1 = getCurrentCharacter();
			if (!isdigit(char1))
				parseError("Failed to get digit for ring closure");
			position++;
			auto char2 = getCurrentCharacter();
			if (!isdigit(char2))
				parseError("Failed to get digit for ring closure");
			position++;
			auto str = string(1, char1) + string(1, char2);
			ringClosureNo = stoi(str);
		}

		auto msg =
				(boost::format("parsed ring closure %d") % ringClosureNo).str();
		REPORT(Reporter::DEBUG) << parseMessage(msg);
		auto keyIter = ringOpenings.find(ringClosureNo);
		if (keyIter != ringOpenings.end()) {
			// ring closing- identify bond type
			auto & ringOpening = keyIter->second;
			assert(ringOpening->getClosureNo() == ringClosureNo);
			if (currentBond != BondType::BondTypeId::UNK
					|| ringOpening->getBondType()
							!= BondType::BondTypeId::UNK) {
				if (currentBond != BondType::BondTypeId::UNK
						&& ringOpening->getBondType()
								!= BondType::BondTypeId::UNK
						&& currentBond != ringOpening->getBondType()) {
					parseError(
							"Ring closure bond type does not match ring opening bond type");
				}
				auto bondType =
						currentBond != BondType::BondTypeId::UNK ?
								currentBond : ringOpening->getBondType();
				addBond(bondType, ringOpening->getAtom(), currentAtom);
				currentBond = BondType::BondTypeId::UNK;
			} else {
				// ring closing default bond
				const auto aromatic = isAromatic(ringOpening->getAtom())
						&& isAromatic(currentAtom);
				const auto & bondType =
						aromatic ?
								BondType::BondTypeId::AR :
								BondType::BondTypeId::SINGLE;
				addBond(bondType, ringOpening->getAtom(), currentAtom);
			}
			ringOpenings.erase(keyIter);
			REPORT(Reporter::DEBUG) << parseMessage("closed ring");
		} else {
			// ring opening
			ringOpenings[ringClosureNo] = make_unique<RingClosureInformation>(
					ringClosureNo, currentAtom, currentBond);
			REPORT(Reporter::DEBUG) << parseMessage("Opened ring");
			currentBond = BondType::BondTypeId::UNK;
		}

		return true;
	}
	return false;
}

void SmilesParser::parseBond() {
	auto currentChar = getCurrentCharacter();

	if (currentChar == '/' || currentChar == '\\') {
		REPORT(Reporter::WARN) << parseMessage("Unsupported stereo bond");
		position++;
		currentChar = getCurrentCharacter();
	}
	auto bondType = BondType::BondTypeId::UNK;
	switch (currentChar) {
	case '-':
		bondType = BondType::BondTypeId::SINGLE;
		position++;
		break;
	case '=':
		bondType = BondType::BondTypeId::DOUBLE;
		position++;
		break;
	case '#':
		bondType = BondType::BondTypeId::TRIPLE;
		position++;
		break;
	case ':':
		bondType = BondType::BondTypeId::AR;
		position++;
		break;
	default:
		break;
	}

	REPORT(Reporter::DEBUG)
			<< parseMessage(
					"parsed bond of type "
							+ BondType::typeFromTypeId(bondType).getName());
	currentBond = bondType;
}

void SmilesParser::parseAtom() {
	auto previousAtom = currentAtom;
	if (getCurrentCharacter() == '[') {
		// read non-organic subset atom type or atom with attributes
		position++;
		// ignore isotropic specifications
		if (isdigit(getCurrentCharacter())) {
			REPORT(Reporter::WARN)
					<< "Ignoring isotropic specification in smiles " + smiles
					<< " at position " << position;
			while (isdigit(getCurrentCharacter())) {
				position++;
			}
		}
		const AtomType & type = parseAtomType(false);
		currentAtom = addAtom(type);
		while (getCurrentCharacter() != ']') {
			parseAtomAttributes();
		}
		assert(getCurrentCharacter() == ']');
		position++;
	} else {
		// atom is single character in organic subset
		const AtomType & type = parseAtomType(true);
		currentAtom = addAtom(type);
	}
	REPORT(Reporter::DEBUG)
			<< (boost::format("Added atom of type %s atom no %d aromatic %s")
					% currentAtom->getAtomType().getName()
					% currentAtom->getAtomNo() % currentAtomAromatic).str();

	if (currentAtomAromatic) {
		aromaticAtomNumbers.insert(currentAtom->getAtomNo());
	}
	if (firstAtomRead && previousAtom != nullptr) {
		// add bond to previous atom
		if (currentBond == BondType::BondTypeId::UNK) {
			currentBond =
					isAromatic(previousAtom) && isAromatic(currentAtom) ?
							BondType::BondTypeId::AR :
							BondType::BondTypeId::SINGLE;
		}
		addBond(currentBond, previousAtom, currentAtom);
	}

	firstAtomRead = true;
}

const bool SmilesParser::isAromatic(const Atom * atom) const {
	return aromaticAtomNumbers.find(atom->getAtomNo())
			!= aromaticAtomNumbers.end();
}

void SmilesParser::parseAtomAttributes() {
	// parse attributes in smiles specification
	auto currentChar = getCurrentCharacter();
	if (currentChar == '+' || currentChar == '-') {
		// charge
		position++;
		auto charge = getCharge(currentChar);
		currentAtom->setFormalCharge(charge);
		currentAtom->setFormalChargeSet(true);
		REPORT(Reporter::DEBUG)
				<< parseMessage(
						(boost::format("added formal charge attribute of %d")
								% charge).str());
		return;
	} else if (currentChar == 'H') {
		// hydrogen count
		int noHydrogens = 1;
		position++;
		auto nextChar = getCurrentCharacter();
		if (isdigit(nextChar)) {
			noHydrogens = nextChar - '0';
			position++;
		}
		addHydrogens(noHydrogens);
		return;
	} else if (currentChar == '@') {
		// stereochemistry - ignore
		position++;
		auto nextChar = getCurrentCharacter();
		if (nextChar == '@') {
			position++;
		}
		REPORT(Reporter::WARN)
				<< parseMessage("Ignoring unsupported stereochemistry");
		return;
	}
	parseError(
			(boost::format("Unknown atom attribute: %s") % currentChar).str());
}

void SmilesParser::addHydrogens(int noHydrogens) {
	addHydrogens(currentAtom, noHydrogens, false);
	REPORT(Reporter::DEBUG)
			<< parseMessage(
					(boost::format("added %d  explicit hydrogens") % noHydrogens).str());
}

void SmilesParser::addHydrogens(const Atom * atom, const int noHydrogens,
		const bool implicit) {
	for (auto i = 0; i < noHydrogens; i++) {
		auto hydrogen = addAtom(
				AtomType::typeFromTypeId(AtomType::AtomTypeId::H), implicit);
		addBond(BondType::BondTypeId::SINGLE, atom, hydrogen);
	}
}

Atom * SmilesParser::addAtom(const AtomType & type, const bool implicit) {
	auto msg = (boost::format("Adding atom %d type %s ") % atomNo
			% type.getName()).str();
	REPORT(Reporter::DEBUG) << parseMessage(msg);
	auto atom = make_unique<Atom>(atomNo, type, type.getName(), implicit);
	auto rtnPtr = atom.get();
	atoms.push_back(move(atom));
	assert(rtnPtr == atoms.at(atoms.size()-1).get());
	atomNo++;
	return rtnPtr;
}

Bond * SmilesParser::addBond(const BondType::BondTypeId typeId,
		const Atom * atom1, const Atom * atom2) {
	assert(atom1 != nullptr);
	assert(atom2 != nullptr);

	const auto & bondType = BondType::typeFromTypeId(typeId);
	auto msg = (boost::format("Adding bond %d type %s between atoms %d and %d")
			% bondNo % bondType.getName() % atom1->getAtomNo()
			% atom2->getAtomNo()).str();
	REPORT(Reporter::DEBUG) << parseMessage(msg);
	auto bond = make_unique<Bond>(bondNo, bondType, atom1, atom2);
	auto bondPtr = bond.get();
	bonds.push_back(move(bond));
	bondNo++;
	return bondPtr;
}

const AtomType & aromaticType(const AtomType & type) {
	// aromatic types for elemental types
	switch (type.getType()) {
	case AtomType::AtomTypeId::C:
		return AtomType::typeFromTypeId(AtomType::AtomTypeId::CAR);
		break;
	case AtomType::AtomTypeId::N:
		return AtomType::typeFromTypeId(AtomType::AtomTypeId::NAR);
		break;
	case AtomType::AtomTypeId::OANY:
		return AtomType::typeFromTypeId(AtomType::AtomTypeId::OAR);
		break;
	default:
		return type;
	}
}

const AtomType & SmilesParser::parseAtomType(const bool organicSubset) {
	char firstCharacter = getCurrentCharacter();
	assert(isalpha(firstCharacter));
	currentAtomAromatic = islower(firstCharacter);
	firstCharacter = toupper(firstCharacter);
	position++;
	string atomSymbol;
	char secondCharacter = getCurrentCharacter();
	if (islower(secondCharacter)) {
		atomSymbol = string(1, firstCharacter) + string(1, secondCharacter);
		if (organicSubset && atomSymbol != "Br" && atomSymbol != "Cl") {
			atomSymbol = string(1, firstCharacter);
		} else {
			position++;
		}
	} else {
		atomSymbol = string(1, firstCharacter);
	}

	const AtomType & elementalType = AtomType::typeFromName(atomSymbol);
	const AtomType & type =
			currentAtomAromatic ? aromaticType(elementalType) : elementalType;

	if (type.getType() == AtomType::AtomTypeId::DU) {
		parseError(
				"Unable to convert atom symbol " + atomSymbol
						+ " to atom type");
	}

	REPORT(Reporter::DEBUG)
			<< parseMessage("parsed atom type " + type.getName());
	return type;
}

void SmilesParser::parseError(const string & message) {
	auto smi = smiles;
	smi.insert(position, "^");
	auto msg = (boost::format("Error parsing smiles %s at position %3d : %s")
			% smi % position % message).str();
	throw runtime_error(msg);
}

const string SmilesParser::parseMessage(const string & message) {
	auto smi = smiles;
	smi.insert(position, "^");
	return (boost::format("Parsing smiles %s at position %3d : %s") % smi
			% position % message).str();

}

char SmilesParser::getCharacter(const size_t pos) const {
	if (pos >= smiles.length()) {
		REPORT(Reporter::TRACE) << "Reading end of smarts";
		return '\0';
	}
	auto ch = smiles.c_str()[pos];
	REPORT(Reporter::TRACE)
			<< boost::format("Read character %s at position %d") % ch % pos;
	return ch;
}

char SmilesParser::getCurrentCharacter() const {
	return getCharacter(position);
}

const int SmilesParser::getCharge(const char chargeChar) {
	assert(chargeChar == '-' || chargeChar == '+');
	auto charge = getInteger();
	if (charge == -1) {
		auto count = countRepeatingCharacter(chargeChar);
		if (count > 0)
			charge = count + 1;
		else
			charge = 1;
	}
	if (chargeChar == '-')
		charge *= -1;
	return charge;
}

const int SmilesParser::getInteger() {
	string str = "";
	char currentChar;
	while (isdigit(currentChar = getCurrentCharacter())) {
		position++;
		str += currentChar;
	}
	if (str == "") {
		return -1;
	}
	return stoi(str);
}

const int SmilesParser::countRepeatingCharacter(const char testChar) {
	int count = 0;
	while (getCurrentCharacter() == testChar) {
		position++;
		count++;
	}
	return count;
}
} /* namespace GarethMol */
