/*
 * SmilesParser.cpp
 *
 *  Created on: Sep 22, 2015
 *      Author: gjones
 */

#include "SmartsParser.h"

#include <boost/format.hpp>
#include <cctype>
#include <iostream>
#include <map>
#include <memory>
#include <vector>

#include "../util/Reporter.h"
#include "Atom.h"
#include "AtomType.h"
#include "BondType.h"
#include "SmartsAstNode.h"

namespace GarethMol {

using namespace std;
using namespace GarethUtil;

SmartsParser::SmartsParser() {
}

SmartsParser::~SmartsParser() {
}

void SmartsParser::parseAtom() {
	auto previousAtom = currentAtom;
	QueryAtomSmartsPtr atomAst = nullptr;
	auto negate = false, firstTerm = true;
	// unless an explicit hydrogen appears first in the atom definition we're matching on a heavy atom
	matchingHeavyAtom = true;
	// check to see if we have a full definition
	if (getCurrentCharacter() == '[') {
		position++;
		vector<unique_ptr<SmartsAstToken<Atom>>> terms;
		while (true) {
			if (getCurrentCharacter() == '!') {
				negate = true;
				position++;
			}
			atomAst = parseAtomPrimitive(firstTerm, negate);
			auto linkage = SmartsParser::getLinkage();
			auto token = make_unique<SmartsAstToken<Atom>>(linkage, atomAst);
			terms.push_back(move(token));
			firstTerm = false;

			if (getCurrentCharacter() == ']') {
				position++;
				break;
			}
		}

		atomAst = SmartsAstToken<Atom>::buildAstTree(terms);
		assert(terms.size() == 1);

	} else {
		// or if it's just a single characater from the organic subset
		atomAst = parseOrganicSubset();
	}

	currentAtom = addAtom(atomAst);

	if (firstAtomRead) {
		addBond(previousAtom, currentAtom);
	}

	firstAtomRead = true;
}

/**
 * Function to match default bond, which is single or aromatic
 *
 * @param bond
 * @return
 */
bool defaultBondMatcher(const Bond & bond) {
	return bond.getBondTypeId() == BondType::BondTypeId::SINGLE
			|| bond.getBondTypeId() == BondType::BondTypeId::AR
			|| bond.getBondTypeId() == BondType::BondTypeId::AM;
}

void SmartsParser::parseBond() {
	// parse bond definion
	vector<unique_ptr<SmartsAstToken<Bond>>> terms;
	while (true) {
		auto negate = false;
		auto currentCharacter = getCurrentCharacter();
		if (currentCharacter == '!') {
			negate = true;
			position++;
			currentCharacter = getCurrentCharacter();
		}

		auto bondAst = parseBondType(negate);
		if (bondAst == nullptr) {
			if (negate) {
				parseError("Negate character with no bond definition");
			}

			if (terms.size() ==0) {
				// no bond definition bound at the current position
				currentBondAst = move(nullptr);
				return;
			}
			break;
		}

		auto linkage = SmartsParser::getLinkage();
		auto token = make_unique<SmartsAstToken<Bond>>(linkage, bondAst);
		terms.push_back(move(token));
	}

	auto bondAst = SmartsAstToken<Bond>::buildAstTree(terms);
	assert(terms.size() == 1);
	// save current bond definition
	currentBondAst = move(bondAst);

}

QueryBondSmartsPtr SmartsParser::parseBondType(const bool negate) {
	auto currentCharacter = getCurrentCharacter();
	if (currentCharacter == '/' || currentCharacter == '\\') {
		parseError("bond stereo chemistry matches are not supported");
	}

	position++;
	auto bondType = BondType::BondTypeId::UNK;
	bool ring = false;
	switch (currentCharacter) {
	case '-':
		bondType = BondType::BondTypeId::SINGLE;
		break;
	case '=':
		bondType = BondType::BondTypeId::DOUBLE;
		break;
	case '#':
		bondType = BondType::BondTypeId::TRIPLE;
		break;
	case ':':
		bondType = BondType::BondTypeId::AR;
		break;
	case '~':
		bondType = BondType::BondTypeId::ANY;
		break;
	case '@':
		bondType = BondType::BondTypeId::ANY;
		ring = true;
		break;
	default:
		position--;
		return nullptr;
		break;
	}

	auto matchFunction = [bondType, ring] (const Bond & bond) {
		if (ring && !bond.isInRing()) {
			return false;
		}
		return BondType::matchType(bond.getBondTypeId(), bondType);
	};

	string description = "Matching bond type "
			+ BondType::typeFromTypeId(bondType).getName()
			+ (ring ? " in ring " : "");
	return make_unique<SmartsAstLeaf<Bond>>(matchFunction, negate, description);
}

const SmartsAstTokenLinkage SmartsParser::getLinkage() {
	char currentCharacter = getCurrentCharacter();
	position++;
	switch (currentCharacter) {
	case '&':
		return SmartsAstTokenLinkage::HIGH_AND;
	case ',':
		return SmartsAstTokenLinkage::OR;
	case ';':
		return SmartsAstTokenLinkage::LOW_AND;
	default:
		position--;
		return SmartsAstTokenLinkage::HIGH_AND;
	}
}

QueryAtomSmartsPtr SmartsParser::parseAtomPrimitive(const bool firstTerm,
		const bool negate) {
	auto currentChar = getCurrentCharacter();
	position++;

	// handle atom primitive terms from the smarts specification'

	if (currentChar == 'D') {
		auto degree = getInteger();
		if (degree == -1) {
			degree = 1;
		}
		const auto matchFunction =
				[degree] (const Atom & atom) {
			return atom.getNeighbourhood().countDegree(false)==degree;
				};
		string description = "Atom of degree " + to_string(degree);
		return make_unique<SmartsAstLeaf<Atom>>(matchFunction, negate,
				description);
	}

	// note H can be a hydrogen count for a heavy atom or an explicit hydrogen atom.
	// As best as I can tell if the H appears immediately in the atom terms then it
	// is a Hydrogen atom, otherwise its a hydrogen count

	if (currentChar == 'H' && !firstTerm) {
		auto nHydrogens = getInteger();
		if (nHydrogens == -1) {
			nHydrogens = 1;
		}
		const auto matchFunction =
				[nHydrogens] (const Atom & atom) {
					return atom.getNeighbourhood().getNoHydrogens()==nHydrogens;
				};
		string description = "Atom has " + to_string(nHydrogens) + " hydrogens";
		return make_unique<SmartsAstLeaf<Atom>>(matchFunction, negate,
				description);
	}         

	else if (currentChar == 'h') {
		auto nHydrogens = getInteger();
		if (nHydrogens == -1) {
			nHydrogens = 1;
		}
		const auto matchFunction =
				[nHydrogens] (const Atom & atom) {
					return atom.getNeighbourhood().countImplicitHydrogens()==nHydrogens;
				};
		string description = "Atom has " + to_string(nHydrogens)
				+ " implicit hydrogens";
		return make_unique<SmartsAstLeaf<Atom>>(matchFunction, negate,
				description);
	}

	if (currentChar == 'R') {
		auto nRings = getInteger();
		if (nRings == -1) {
			const auto matchFunction =
					[nRings] (const Atom & atom) {
						return atom.getNeighbourhood().countRings()>0;
					};
			string description = "Atom is in in any ring";
			return make_unique<SmartsAstLeaf<Atom>>(matchFunction, negate,
					description);
		}
		const auto matchFunction =
				[nRings] (const Atom & atom) {
					return atom.getNeighbourhood().countRings()==nRings;
				};
		string description = "Atom is in  " + to_string(nRings) + " SSSR rings";
		return make_unique<SmartsAstLeaf<Atom>>(matchFunction, negate,
				description);
	}

	if (currentChar == 'r') {
		auto ringSize = getInteger();
		if (ringSize == -1) {
			const auto matchFunction =
					[ringSize] (const Atom & atom) {
						return atom.getNeighbourhood().countRings()>0;
					};
			string description = "Atom is in a ring";
			return make_unique<SmartsAstLeaf<Atom>>(matchFunction, negate,
					description);
		}
		const auto matchFunction =
				[ringSize] (const Atom & atom) {
					return atom.getNeighbourhood().inRingOfSize(ringSize);
				};
		string description = "Atom is in a ring of size " + to_string(ringSize);
		return make_unique<SmartsAstLeaf<Atom>>(matchFunction, negate,
				description);
	}

	else if (currentChar == 'v') {
		auto valence = getInteger();
		if (valence == -1) {
			valence = 1;
		}
		const auto matchFunction =
				[valence] (const Atom & atom) {
					return atom.getNeighbourhood().totalBondOrder()==valence;
				};
		string description = "Atom has valence " + to_string(valence);
		return make_unique<SmartsAstLeaf<Atom>>(matchFunction, negate,
				description);
	}

	else if (currentChar == 'X') {
		auto nConnections = getInteger();
		if (nConnections == -1) {
			nConnections = 1;
		}
		const auto matchFunction =
				[nConnections] (const Atom & atom) {
					return atom.getNeighbourhood().getAtomNeighbours().size()==static_cast<size_t>(nConnections);
				};
		string description = "Atom has " + to_string(nConnections)
				+ " connections";
		return make_unique<SmartsAstLeaf<Atom>>(matchFunction, negate,
				description);
	}

	else if (currentChar == 'x') {
		auto nRingConnections = getInteger();
		if (nRingConnections == -1) {
			nRingConnections = 1;
		}
		const auto matchFunction =
				[nRingConnections] (const Atom & atom) {
					return atom.getNeighbourhood().countRingConnections()==nRingConnections;
				};
		string description = "Atom has " + to_string(nRingConnections)
				+ " ring connections";
		return make_unique<SmartsAstLeaf<Atom>>(matchFunction, negate,
				description);
	}

	else if (currentChar == '+' || currentChar == '-') {
		return move(getCharge(negate, currentChar));
	}

	else if (currentChar == '@') {
		parseError("Chirality matches are not supported");
	}

	else if (isdigit(currentChar)) {
		parseError("Atomic weight matches are not supported");
	}

	else if (currentChar == '$') {
		string recursiveSmarts = extractRecursiveSmarts();
		SmartsParser recursiveSmartsParser;
		auto recursiveSmartsQueryMolecule = recursiveSmartsParser.parseSmarts(recursiveSmarts);
		return make_unique<SmartsAstRecursiveLeaf>(recursiveSmartsQueryMolecule, negate);
		parseError("Recursive smarts are not supported");
	}

	position--;

	// not a primitive, but an atom type definition- an initial H will fall through to here
	return parseAtomType(firstTerm, negate);

}

string SmartsParser::extractRecursiveSmarts() {
	auto currentChar = getCurrentCharacter();
	assert(currentChar == '(');
	auto nOpen = 1;
	position++;
	auto start = position;
	while(true) {
		currentChar = getCurrentCharacter();
		position++;
		if (currentChar =='(') {
			nOpen++;
		}
		else if (currentChar == ')') {
			nOpen--;
		}
		else if (currentChar == '\0') {
			break;
		}
		if (nOpen == 0) {
			auto recursiveSmarts = smarts.substr(start, position-start-1);
			return recursiveSmarts;
		}
	}
	parseError("Unable to find termination of recursive smarts");
	return string("");
}

QueryAtomSmartsPtr SmartsParser::parseOrganicSubset() {
	// atoms from the smiles orginic subset and aromatic/aliphatic(a/A)
	auto currentChar = getCurrentCharacter();
	position++;
	auto nextChar = getCurrentCharacter();
	bool aromatic = islower(currentChar);
	currentChar = toupper(currentChar);
	AtomType::AtomTypeId atomType = AtomType::AtomTypeId::ATM_NONE;
	if (currentChar == 'C' && nextChar == 'l') {
		position++;
		atomType = AtomType::AtomTypeId::CL;
	} else if (currentChar == 'B' && nextChar == 'r') {
		position++;
		atomType = AtomType::AtomTypeId::BR;
	} else if (currentChar == 'C') {
		atomType = AtomType::AtomTypeId::C;
	} else if (currentChar == 'B') {
		atomType = AtomType::AtomTypeId::B;
	} else if (currentChar == 'N') {
		atomType = AtomType::AtomTypeId::N;
	} else if (currentChar == 'O') {
		atomType = AtomType::AtomTypeId::OANY;
	} else if (currentChar == 'P') {
		atomType = AtomType::AtomTypeId::P;
	} else if (currentChar == 'S') {
		atomType = AtomType::AtomTypeId::S;
	} else if (currentChar == 'F') {
		atomType = AtomType::AtomTypeId::F;
	} else if (currentChar == 'I') {
		atomType = AtomType::AtomTypeId::I;
	} else if (currentChar == '*') {
		atomType = AtomType::AtomTypeId::ANY;
		const auto matchFunction = [] (const Atom & atom) {
			return true;
		};
		string description = "Atom is any type";
		return make_unique<SmartsAstLeaf<Atom>>(matchFunction, false,
				description);
	} else if (currentChar == 'A' && aromatic) {
		const auto matchFunction = [] (const Atom & atom) {
			return atom.getNeighbourhood().isAromatic();
		};
		string description = "Atom is any aromatic type";
		return make_unique<SmartsAstLeaf<Atom>>(matchFunction, false,
				description);
	} else if (currentChar == 'A') {
		const auto matchFunction = [] (const Atom & atom) {
			return !atom.getNeighbourhood().isAromatic();
		};
		string description = "Atom is any aliphatic type";
		return make_unique<SmartsAstLeaf<Atom>>(matchFunction, false,
				description);
	}

	if (atomType != AtomType::AtomTypeId::ATM_NONE) {
		const auto matchFunction = [atomType, aromatic] (const Atom & atom) {
			if (! AtomType::matchType(atomType, atom.getAtomTypeId())) {
				return false;
			}
			return aromatic == atom.getNeighbourhood().isAromatic();
		};
		string description = "Atom is " + AtomType::nameFromTypeId(atomType)
				+ " and " + (aromatic ? "aromatic" : "aliphatic");
		return make_unique<SmartsAstLeaf<Atom>>(matchFunction, false,
				description);
	}

	parseError("Failed to parse atom definition");
	return nullptr;
}

QueryAtomSmartsPtr SmartsParser::parseAtomType(const bool firstTerm,
		const bool negate) {

	// more general atom definition for types in the full atom definition

	char currentChar = getCurrentCharacter();
	position++;

	bool aromatic = false;
	bool aliphatic = false;
	AtomType::AtomTypeId atomType = AtomType::AtomTypeId::ATM_NONE;

	string atomSymbol = "";
	if (currentChar == '*') {
		atomType = AtomType::AtomTypeId::ANY;
	} else if (currentChar == 'a') {
		atomType = AtomType::AtomTypeId::ANY;
		aromatic = true;
	} else if (currentChar == 'A') {
		atomType = AtomType::AtomTypeId::ANY;
		aliphatic = true;
	} else if (islower(currentChar)) {
		aromatic = true;
		atomSymbol = string(1, toupper(currentChar));
	} else if (isupper(currentChar)) {
		char secondCharacter = getCurrentCharacter();
		// two tests here - if the second character is a v then we are looking at
		// the bond order primitive not an atom type.  If the second character
		// is r that's a ring check, unless it's in Br
		if (islower(secondCharacter) && secondCharacter != 'v'
				&& (currentChar == 'B' || secondCharacter != 'r')) {
			position++;
			atomSymbol = string(1, currentChar) + string(1, secondCharacter);
		} else {
			if (currentChar == 'N' || currentChar == 'C' || currentChar == 'S'
					|| currentChar == 'O')
				aliphatic = true;
			atomSymbol = string(1, currentChar);
		}
	} else if (currentChar == '#') {
		const int atomicNumber = getInteger();
		const auto & type = AtomType::typeFromAtomicNumber(atomicNumber);
		atomType = type.getType();
	}

	// get atom for symbol
	if (atomSymbol != "") {
		const AtomType & type = AtomType::typeFromName(atomSymbol);
		if (type.getType() == AtomType::AtomTypeId::DU) {
			parseError(
					"Unable to convert atom symbol " + atomSymbol
							+ " to atom type");
		}
		atomType = type.getType();
	}

	// check for explicit hydrogen match
	if (atomType == AtomType::AtomTypeId::H) {
		matchingHeavyAtom = false;
	}

	// create AST leaf node for atom
	if (atomType != AtomType::AtomTypeId::ATM_NONE || aromatic || aliphatic) {
		const auto matchFunction =
				[atomType, aromatic, aliphatic] (const Atom & atom) {
					if (atomType != AtomType::AtomTypeId::ATM_NONE) {
						if (! AtomType::matchType(atomType, atom.getAtomTypeId()))
						return false;
					}
					if (aromatic && !atom.getNeighbourhood().isAromatic()) {
						return false;
					}
					if (aliphatic && atom.getNeighbourhood().isAromatic()) {
						return false;
					}
					return true;
				};
		string description = "Atom is " + AtomType::nameFromTypeId(atomType)
				+ (aromatic ? " and aromatic " : "")
				+ (aliphatic ? " and aliphatic " : "");
		return make_unique<SmartsAstLeaf<Atom>>(matchFunction, negate,
				description);
	}

	parseError("Failed to parse atom definition");
	return nullptr;
}

QueryAtomSmartsPtr SmartsParser::getCharge(const bool negate,
		const char chargeChar) {
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
	const auto matchFunction = [charge] (const Atom & atom) {
		return atom.getFormalCharge()==charge;
	};
	string description = "Atom has formal charge " + to_string(charge);
	return make_unique<SmartsAstLeaf<Atom>>(matchFunction, negate, description);
}

const int SmartsParser::getInteger() {
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

const int SmartsParser::countRepeatingCharacter(const char testChar) {
	int count = 0;
	while (getCurrentCharacter() == testChar) {
		position++;
		count++;
	}
	return count;
}

unique_ptr<QueryMolecule> SmartsParser::parseSmarts(const string & smarts_) {

	// public function to create query molecule

	// identify any name for the smarts
	auto pos = smarts_.find(' ');
	string name = "Smarts";
	smarts = pos == string::npos ? smarts_ : smarts_.substr(0, pos);
	if (pos != string::npos) {
		name = smarts_.substr(pos + 1, string::npos);
	}

	// initialize
	position = 0;
	atomNo = 0;
	bondNo = 0;
	atoms.clear();
	bonds.clear();
	ringOpenings= {};
	openBranches = {};
	atoms.reserve(smarts.length() * 2);
	bonds.reserve(smarts.length() * 2);
	aromaticAtomNumbers = {};
	currentAtom = nullptr;
	currentBondAst = nullptr;
	firstAtomRead = false;

	// parse
	while (position < smarts.length()) {
		parseAtomEntry();
	}
	auto molecule = make_unique<QueryMolecule>(name + " " + smarts, atoms,
			bonds);

	return molecule;
}

void SmartsParser::parseAtomEntry() {
	// parses atom and related bond entries

	if (firstAtomRead) {
		// branches and bond to previous atom
		if (getCurrentCharacter() == '.')
			parseError("Disconnected structures are not supported");
		auto openBranch = parseOpenBranch();
		if (openBranch && currentBondAst != nullptr) {
			parseError("Unused bond definition at open branch");
		}
		if (currentBondAst == nullptr)
			parseBond();
	}
	// get atom entry
	parseAtom();
	// handle ring closures/openings
	auto readingRingClosure = true;
	while (readingRingClosure) {
		// we don't know (yet) whether or not this bond is for a ring closure or for
		// the bond to the next atom
		parseBond();
		readingRingClosure = parseRingClosure();
	}
	// any close branch
	parseCloseBranch();
}

void SmartsParser::parseCloseBranch() {
	while (getCurrentCharacter() == ')') {
		position++;
		if (openBranches.size() == 0) {
			parseError("No matching open branch for close");
		}
		REPORT(Reporter::TRACE) << parseMessage("parsed close branch");
		currentAtom = openBranches.back();
		openBranches.pop_back();
	}
}

bool SmartsParser::parseOpenBranch() {
	if (getCurrentCharacter() == '(') {
		position++;
		openBranches.push_back(currentAtom);
		REPORT(Reporter::TRACE) << parseMessage("parsed open branch");
		return true;
	}
	return false;
}

bool SmartsParser::parseRingClosure() {
	auto currentChar = getCurrentCharacter();
	if (isdigit(currentChar) || currentChar == '%') {
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
		REPORT(Reporter::TRACE) << parseMessage(msg);
		auto keyIter = ringOpenings.find(ringClosureNo);
		if (keyIter != ringOpenings.end()) {
			auto ringOpeningAtom = keyIter->second;
			addBond(ringOpeningAtom, currentAtom);
			ringOpenings.erase(keyIter);
			REPORT(Reporter::TRACE) << parseMessage("closed ring");
		} else {
			ringOpenings[ringClosureNo] = currentAtom;
			REPORT(Reporter::TRACE) << parseMessage("Opened ring");
		}

		return true;
	}
	return false;
}

QueryAtom * SmartsParser::addAtom(QueryAtomSmartsPtr & atomSmarts) {
	auto atom = make_unique<QueryAtom>(atomNo, matchingHeavyAtom, atomSmarts);
	auto msg = (boost::format("Adding query atom %d ") % atomNo).str();
	REPORT(Reporter::TRACE) << parseMessage(msg);
	atomNo++;
	auto atomPtr = atom.get();
	atoms.push_back(move(atom));
	return atomPtr;
}

QueryBond * SmartsParser::addBond(const QueryAtom * atom1,
		const QueryAtom * atom2) {
	assert(atom1 != nullptr);
	assert(atom2 != nullptr);

	if (currentBondAst == nullptr) {
		currentBondAst = make_unique<SmartsAstLeaf<Bond>>(defaultBondMatcher,
				false, "Default bond is single or aromatic");
	}
	auto bond = make_unique<QueryBond>(bondNo, currentBondAst, atom1, atom2);
	auto msg = (boost::format("Adding bond %d between atoms %d and %d") % bondNo
			% atom1->getAtomNo() % atom2->getAtomNo()).str();
	REPORT(Reporter::TRACE) << parseMessage(msg);
	bondNo++;
	currentBondAst = nullptr;
	auto bondPtr = bond.get();
	bonds.push_back(move(bond));
	return bondPtr;
}

void SmartsParser::parseError(const string & message) {
	auto smi = smarts;
	smi.insert(position, "^");
	auto msg = (boost::format("Error parsing smarts %s at position %3d : %s")
			% smi % position % message).str();
	throw runtime_error(msg);
}

const string SmartsParser::parseMessage(const string & message) {
	auto smi = smarts;
	smi.insert(position, "^");
	return (boost::format("Parsing smarts %s at position %3d : %s") % smi
			% position % message).str();

}

char SmartsParser::getCharacter(const size_t pos) const {
	if (pos >= smarts.length()) {
		REPORT(Reporter::TRACE) << "Reading end of smarts";
		return '\0';
	}
	auto ch = smarts.c_str()[pos];
	REPORT(Reporter::TRACE)
			<< boost::format("Read character %s at position %d") % ch % pos;
	return ch;
}

char SmartsParser::getCurrentCharacter() const {
	return getCharacter(position);
}

} /* namespace GarethMol */
