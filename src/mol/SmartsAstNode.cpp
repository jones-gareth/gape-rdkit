/*
 * SmartsAstNode.cpp
 *
 *  Created on: Jun 2, 2016
 *      Author: gareth
 */

#include "QueryMolecule.h"
#include "SmartsAstNode.h"
#include "SubstructureSearch.h"

namespace GarethMol {

using namespace std;

SmartsAstRecursiveLeaf::SmartsAstRecursiveLeaf(
		unique_ptr<QueryMolecule> & recursiveSmartsPattern_, const bool negate_) :
		SmartsAstNode<Atom>(SmartsAstNodeType::LEAF), recursiveSmartsPattern(
				move(recursiveSmartsPattern_)), negate(negate_) {
}

SmartsAstRecursiveLeaf::~SmartsAstRecursiveLeaf()  {
	;
}
/**
 * Match an atom against a recursive smarts
 *
 * @param target
 * @return
 */
bool SmartsAstRecursiveLeaf::evaluate(const Atom & target,
		const Molecule & targetMolecule) const {
	subSearch->fixMapping(0, target.getAtomNo());
	auto nIsomorphisms = subSearch->compare();
	auto value = nIsomorphisms > 0;
	subSearch->clearFixedMappings();
	auto rtn = negate ? !value : value;
	REPORT(Reporter::TRACE) << boolalpha << "Evaluated AST "
			<< typeid(Atom).name() << " leaf node:"
			<< (negate ? " [negate value] " : " ") << rtn << " '" << description
			<< "'" << " against target " << target.info();
	return rtn;
}

void SmartsAstRecursiveLeaf::initialize(const Molecule & targetMolecule) {
	subSearch = make_unique<SubstructureSearch>(*recursiveSmartsPattern,
			targetMolecule);
}

}

