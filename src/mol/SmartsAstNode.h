/*
 * SmartsAst.h
 *
 *  Created on: Sep 26, 2015
 *      Author: gjones
 *
 * Classes to support the representation of Abstract Syntax Trees (ASTs) for query
 * atoms and bonds
 */

#ifndef SRC_MOL_SMARTSASTNODE_H_
#define SRC_MOL_SMARTSASTNODE_H_

#include "../util/Reporter.h"
#include <functional>
#include <type_traits>
#include "Atom.h"
#include "Bond.h"

namespace GarethMol {

using namespace std;
using namespace GarethUtil;
class QueryMolecule;
class SubstructureSearch;

enum class SmartsAstNodeType {
	OR, AND, LEAF
};

/**
 * Base class for a AST node.  Target is an atom or bond
 */
template<class Target>
class SmartsAstNode {
public:
	static_assert(std::is_same<Target, Atom>::value || std::is_same<Target, Bond>::value,
			"Smarts AST node is available for Atom and Bond classes only");

	/**
	 * Construct using child nodes.
	 *
	 * @param nodeType_
	 * @param child1_
	 * @param child2_
	 */
	SmartsAstNode(const SmartsAstNodeType nodeType_,
			unique_ptr<SmartsAstNode<Target>> & child1_,
			unique_ptr<SmartsAstNode<Target>> & child2_) :
			nodeType(nodeType_), child1(move(child1_)), child2(move(child2_)) {
	}

	virtual ~SmartsAstNode() {
	}

	SmartsAstNode(const SmartsAstNode & rhs) = delete;
	SmartsAstNode & operator =(const SmartsAstNode & rhs) = delete;
	SmartsAstNode(SmartsAstNode && rhs) = delete;
	SmartsAstNode & operator =(SmartsAstNode && rhs) = delete;

	/**
	 * Attempts to match the target (atom or bond) against the tree.  Evaluates
	 * down the tree until we test leaf nodes.
	 *
	 * @param target
	 * @return
	 */
	virtual bool evaluate(const Target & target,
			const Molecule & targetMolecule) const {
		bool value = false;
		REPORT(Reporter::TRACE) << "Processing compound AST node type "
				<< nodeTypeToString(nodeType);
		switch (nodeType) {
		case SmartsAstNodeType::AND:
			value = child1->evaluate(target, targetMolecule)
					&& child2->evaluate(target, targetMolecule);
			break;
		case SmartsAstNodeType::OR:
			value = child1->evaluate(target, targetMolecule)
					|| child2->evaluate(target, targetMolecule);
			break;
		default:
			throw runtime_error("Leaf node evaluate not implemented");
		}
		return value;
	}

	virtual void initialize(const Molecule & targetMolecule) {
		switch (nodeType) {
		case SmartsAstNodeType::AND:
		case SmartsAstNodeType::OR:
			child1->initialize(targetMolecule);
			child2->initialize(targetMolecule);
			break;
		default:
			throw runtime_error("Leaf node initialize not implemented");
		}
	}

protected:
	SmartsAstNode(const SmartsAstNodeType nodeType_) :
			nodeType(nodeType_), child1(nullptr), child2(nullptr) {
	}

private:
	const SmartsAstNodeType nodeType;
	unique_ptr<SmartsAstNode<Target>> child1, child2;

	const string nodeTypeToString(SmartsAstNodeType nodeType) const {
		switch (nodeType) {
		case SmartsAstNodeType::AND:
			return "and";
		case SmartsAstNodeType::OR:
			return "or";
		default:
			throw runtime_error("Unknown smarts AST node type");

		}
	}
};

/**
 * Subclass of AST node for a leaf
 */
template<typename Target>
class SmartsAstLeaf: public SmartsAstNode<Target> {
public:
	static_assert(std::is_same<Target, Atom>::value || std::is_same<Target, Bond>::value,
			"Smarts AST leaf is available for Atom and Bond classes only");

	/**
	 * Construct using a function that tests the target (atom or bond).
	 *
	 * @param critiaFunc_ target function
	 * @param negate_ Set true to apply not logic to target function
	 * @param description_
	 */
	explicit SmartsAstLeaf(const function<bool(const Target &)> & critiaFunc_,
			const bool negate_, const string & description_) :
			SmartsAstNode<Target>(SmartsAstNodeType::LEAF), critiaFunc(
					critiaFunc_), negate(negate_), description(description_) {
	}

	SmartsAstLeaf(const SmartsAstLeaf & rhs) = delete;
	SmartsAstLeaf & operator =(const SmartsAstLeaf & rhs) = delete;
	SmartsAstLeaf(SmartsAstLeaf && rhs) = delete;
	SmartsAstLeaf & operator =(SmartsAstLeaf && rhs) = delete;

	/**
	 * Evaluate using the target function
	 *
	 * @param target
	 * @return
	 */
	virtual bool evaluate(const Target & target,
			const Molecule & targetMolecule) const override {
		auto value = critiaFunc(target);
		auto rtn = negate ? !value : value;
		REPORT(Reporter::TRACE) << boolalpha << "Evaluated AST "
				<< typeid(Target).name() << " leaf node:"
				<< (negate ? " [negate value] " : " ") << rtn << " '"
				<< description << "'" << " against target " << target.info();
		return rtn;
	}

	virtual void initialize(const Molecule & targetMolecule) override {
	}

	virtual ~SmartsAstLeaf() override {
	}

private:
	const function<bool(const Target &)> critiaFunc;
	const bool negate;
	const string description;
};

/**
 * Subclass of AST node for a recursive smarts leaf
 */

class SmartsAstRecursiveLeaf: public SmartsAstNode<Atom> {
public:

	/**
	 * Construct using a function that tests the target (atom or bond).
	 *
	 * @param critiaFunc_ target function
	 * @param negate_ Set true to apply not logic to target function
	 * @param description_
	 */
	explicit SmartsAstRecursiveLeaf(
			unique_ptr<QueryMolecule> & recursiveSmartsPattern_,
			const bool negate_);

	SmartsAstRecursiveLeaf(const SmartsAstRecursiveLeaf & rhs) = delete;
	SmartsAstRecursiveLeaf & operator =(const SmartsAstRecursiveLeaf & rhs) = delete;
	SmartsAstRecursiveLeaf(SmartsAstRecursiveLeaf && rhs) = delete;
	SmartsAstRecursiveLeaf & operator =(SmartsAstRecursiveLeaf && rhs) = delete;

	 virtual void initialize(const Molecule & targetMolecule) override;

	virtual bool evaluate(const Atom & target,
			const Molecule & targetMolecule) const override;

	virtual ~SmartsAstRecursiveLeaf() override;

private:
	unique_ptr<QueryMolecule> recursiveSmartsPattern;
	const bool negate;
	const string description = "Appling recursive smarts to atom";
	unique_ptr<SubstructureSearch> subSearch;
};

} /* namespace GarethMol */

#endif /* SRC_MOL_SMARTSASTNODE_H_ */
