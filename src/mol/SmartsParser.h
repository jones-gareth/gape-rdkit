/*
 *
 * Classes to implement a smarts parser.  Does not implement the following features:
 *
 * Isotropic features
 * Bond steochemistry
 * Atom chirality
 * Recursive smarts
 * Disconnected queries/multiple fragments
 *
 * Matching performed using Abstract Syntax Trees for both bonds and atoms.
 *
 * Atom and bond elements in the smarts string are parsed into tokens linked by
 * logical operators which are then converted to ASTs that are representations of
 * query atoms and bonds.
 *
 * See also QueryMolecule.h and SmartsAstNode.h
 *
 * SmartsParser.h
 *
 *  Created on: Sep 22, 2015
 *      Author: gjones
 */

#ifndef SRC_MOL_SMARTSPARSER_H_
#define SRC_MOL_SMARTSPARSER_H_

#include <set>
#include <string>
#include "QueryMolecule.h"

namespace GarethMol {

/**
 * Logical operators for us in tokenization (not including NOT)
 */
enum class SmartsAstTokenLinkage {
	HIGH_AND, OR, LOW_AND
};

/**
 * Template for an AST token. Target is either an Atom or Bond
 */
template<typename Target>
class SmartsAstToken {
public:
	static_assert(std::is_same<Target, Atom>::value || std::is_same<Target, Bond>::value,
			"Smarts AST token is available for Atom and Bond classes only");

	/**
	 * Create a new token with specified linkage and AST node
	 *
	 * @param linkage_
	 * @param smartsAst_
	 */
	SmartsAstToken(SmartsAstTokenLinkage linkage_,
			unique_ptr<SmartsAstNode<Target>> & smartsAst_) :
			linkage(linkage_), smartsAst(move(smartsAst_)) {
	}

	virtual ~SmartsAstToken() {
	}

	SmartsAstToken(const SmartsAstToken & rhs) = delete;
	SmartsAstToken & operator =(const SmartsAstToken & rhs) = delete;
	SmartsAstToken(SmartsAstToken && rhs) = delete;
	SmartsAstToken & operator =(SmartsAstToken && rhs) = delete;

	/**
	 * Given a list of tokens construct an AST tree for node matching
	 *
	 * @param tokens
	 * @return
	 */
	static unique_ptr<SmartsAstNode<Target>> buildAstTree(
			vector<unique_ptr<SmartsAstToken<Target>>>& tokens);

	const SmartsAstTokenLinkage getLinkage() const {
		return linkage;
	}

	unique_ptr<SmartsAstNode<Target> >& getSmartsAst() {
		return smartsAst;
	}

private:
	const SmartsAstTokenLinkage linkage;
	unique_ptr<SmartsAstNode<Target>> smartsAst;
	static bool mergeFirstToken(const SmartsAstTokenLinkage linkage,
	vector<unique_ptr<SmartsAstToken<Target>>> & tokens);

};

/**
 * Implementation class for parsing smarts
 */
class SmartsParser {
public:

	SmartsParser();
	virtual ~SmartsParser();

	SmartsParser(const SmartsParser & rhs) = delete;
	SmartsParser & operator =(const SmartsParser & rhs) = delete;
	SmartsParser(SmartsParser && rhs) = delete;
	SmartsParser & operator =(SmartsParser && rhs) = delete;

	/**
	 * Parse smarts and return Query Molecule
	 *
	 * @param smarts_
	 * @return
	 */
	unique_ptr<QueryMolecule> parseSmarts(const string & smarts_);

private:
	// current position
	size_t position = 0;
	// smarts we're parsing
	string smarts { "" };
	// Query Molecule under construction
	vector<QueryAtomPtr> atoms { };
	vector<QueryBondPtr> bonds { };
	bool firstAtomRead = false;
	QueryAtom * currentAtom = nullptr;
	size_t atomNo = 0, bondNo = 0;
	QueryBondSmartsPtr currentBondAst = nullptr;
	// branches stack
	vector<QueryAtom *> openBranches { };
	// ring openings list
	map<size_t, QueryAtom *> ringOpenings { };
	set<size_t> aromaticAtomNumbers { };
	// normally we only match on heavy atoms, but this is not true if the current
	// atom explicitly matches a Hydrogen
	bool matchingHeavyAtom = true;

	/**
	 * Parse logical operator from smarts
	 * @return
	 */
	const SmartsAstTokenLinkage getLinkage();

	/**
	 * Parse atom type definition
	 *
	 * @param firstTerm
	 * @param negate
	 * @return
	 */
	QueryAtomSmartsPtr parseAtomType(const bool firstTerm, const bool negate);

	/**
	 * Parse simple atom type from organic subset
	 *
	 * @return
	 */
	QueryAtomSmartsPtr parseOrganicSubset();

	/**
	 * Parse general atom type terms and atom attributes (one at a time)
	 *
	 * @param firstTerm
	 * @param negate
	 * @return
	 */
	QueryAtomSmartsPtr parseAtomPrimitive(const bool firstTerm,
			const bool negate);

	/**
	 * Extract the smarts pattern for recursive smarts
	 *
	 * @return
	 */
	string extractRecursiveSmarts();

	/**
	 * Parse a bond definition
	 *
	 * @param negate
	 * @return
	 */
	QueryBondSmartsPtr parseBondType(const bool negate);

	/**
	 * Get next integer from smarts string
	 * @return
	 */
	const int getInteger();

	/**
	 * Get character at position from smarts string
	 * @param pos
	 * @return
	 */
	char getCharacter(const size_t pos) const;

	/**
	 * Get current character from smarts string
	 * @return
	 */
	char getCurrentCharacter() const;

	/**
	 * Handle errors (throws runtime_errror)
	 * @param message
	 */
	void parseError(const string & message);

	/**
	 * Constructs an message, including parse position
	 *
	 * @param message
	 * @return
	 */
	const string parseMessage(const string & message);

	/**
	 * Parses entire atom definition (branches, rings, atom and bond)
	 */
	void parseAtomEntry();

	/**
	 * Handles open branch
	 * @return
	 */
	bool parseOpenBranch();

	/**
	 * Parses bond definition
	 */
	void parseBond();

	/**
	 * Parses atom terms
	 */
	void parseAtom();

	/**
	 * Parses ring closure/opening
	 * @return
	 */
	bool parseRingClosure();

	/**
	 * Handles close branch
	 */
	void parseCloseBranch();

	/**
	 * Given the atom AST adds a atom to the building query molecule.
	 * @param smartsAst
	 * @return
	 */
	QueryAtom * addAtom(QueryAtomSmartsPtr & smartsAst);

	/**
	 * Given two atoms AST (and the current bond AST stored in the class)
	 * adds a bond to the building query molecule.
	 *
	 * @param atom1
	 * @param atom2
	 * @return
	 */
	QueryBond * addBond(const QueryAtom * atom1, const QueryAtom * atom2);

	/**
	 * Gets a count for a character from the smarts string
	 *
	 * @param testChar
	 * @return
	 */
	const int countRepeatingCharacter(const char testChar);

	/**
	 * Gets a charge from the smarts string
	 * @param negate
	 * @param chargeChar
	 * @return
	 */
	QueryAtomSmartsPtr getCharge(const bool negate, const char chargeChar);
};

/**
 * Given a series od tokens contructs an AST tree.
 *
 * @param tokens
 * @return
 */
template<typename Target>
unique_ptr<SmartsAstNode<Target>> SmartsAstToken<Target>::buildAstTree(
		vector<unique_ptr<SmartsAstToken<Target>>>& tokens) {
	// logical operators in order of precidence
	auto operators = vector<SmartsAstTokenLinkage> {
		SmartsAstTokenLinkage::HIGH_AND, SmartsAstTokenLinkage::OR,
		SmartsAstTokenLinkage::LOW_AND};
	// at each level of precidence
	for (auto linkage: operators) {
		// merge two adjacent tokens with linkage at that precedence into one token
		while(mergeFirstToken(linkage, tokens)) {
			;
		}
	}
	// Should have only one token at the end containing the AST
	assert(tokens.size()==1);
	return move(tokens.at(0)->getSmartsAst());
}

/**
 * Look for two tokens with this linkage and join them
 * @param linkage
 * @param tokens
 * @return
 */
template<typename Target>
bool SmartsAstToken<Target>::mergeFirstToken(
		const SmartsAstTokenLinkage linkage,
		vector<unique_ptr<SmartsAstToken<Target>>>&tokens) {
			// find the next token with this linkage
			auto matcher =
			[linkage] (const unique_ptr<SmartsAstToken<Target>> & token) {return token->getLinkage() == linkage;};
			auto iter = find_if(tokens.begin(), tokens.end(), matcher);
			if (iter == tokens.end() || next(iter) == tokens.end())
			return false;

			// identify two tokens to join
			auto & firstToken = *iter;
			auto & nextToken = *next(iter);
			// join into new node
			SmartsAstNodeType nodeType;
			switch (linkage) {
				case SmartsAstTokenLinkage::HIGH_AND:
				case SmartsAstTokenLinkage::LOW_AND:
				nodeType=SmartsAstNodeType::AND;
				break;
				case SmartsAstTokenLinkage::OR:
				nodeType=SmartsAstNodeType::OR;
				break;
			}
			auto newNode = make_unique<SmartsAstNode<Target>>(nodeType,
			firstToken->getSmartsAst(), nextToken->getSmartsAst());
			// replace first token with new token and remove second
			auto newToken = make_unique<SmartsAstToken<Target>>(nextToken->getLinkage(),
			newNode);
			*iter = move(newToken);
			tokens.erase(next(iter));

			return true;
		}

	}
	/* namespace GarethMol */

#endif /* SRC_MOL_SMILESPARSER_H_ */
