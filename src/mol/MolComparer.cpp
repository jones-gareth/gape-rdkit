/*
 * MolComparer.cpp
 *
 *  Created on: Sep 3, 2015
 *      Author: gjones
 */

#include "MolComparer.h"
#include "../util/ListProcessor.h"
#include "../util/Reporter.h"
#include "../util/Util.h"

namespace GarethMol {

using namespace GarethUtil;

void buildAdjacencyMatrix(const Molecule & molecule,
        const std::vector<const Atom *> & atoms,
        Array2D<bool> & adjacencyMatrix,
        std::vector<size_t> & atomConnections) {
    const auto size = atoms.size();
    assert(adjacencyMatrix.getRows() == size);
    assert(adjacencyMatrix.getColumns() == size);
    assert(atomConnections.size() == size);

    for (size_t i = 0; i < size; i++) {
        const auto & atom1 = atoms.at(i);
        for (size_t j = i + 1; j < size; j++) {
            const auto & atom2 = atoms.at(j);
            auto bondOptional = molecule.isBonded(*atom1, *atom2);
            bool isBonded = bondOptional.is_initialized();
            adjacencyMatrix.set(i, j, isBonded);
            adjacencyMatrix.set(j, i, isBonded);

            // don't base atom connection count based on the current matching atoms-
            // even if we are matching heavy atoms only we want to include hydrogens
            // in the total connection count.

//			if (isBonded) {
//				++atomConnections.at(i);
//				++atomConnections.at(j);
//			}

        }
    }

    // for the atom connections count all connections, including hydrogens,
    // even if we are only matching on the heavy atoms
    // (but only include real atoms - skip LP, Du).
    for (size_t i = 0; i < size; i++) {
        atomConnections.at(i) = 0;
        const auto & atom = atoms.at(i);
        if (!AtomType::isRealAtom(atom->getAtomTypeId()))
            continue;
        for (const auto & bond : molecule.getBonds()) {
            const auto & atom1 = bond->getAtom1();
            const auto & atom2 = bond->getAtom2();
            if (atom1.getAtomNo() == atom->getAtomNo()
                    && AtomType::isRealAtom(atom2.getAtomTypeId())) {
                ++atomConnections.at(i);
            } else if (atom2.getAtomNo() == atom->getAtomNo()
                    && AtomType::isRealAtom(atom1.getAtomTypeId())) {
                ++atomConnections.at(i);
            }
        }
    }
}

size_t MolComparer::compare() {

    // find all heavy atoms in the molecule
    queryAtoms = findAtoms(queryMolecule);
    targetAtoms = findAtoms(targetMolecule);

    REPORT(Reporter::TRACE) << "query molecule no of atoms "
            << queryMolecule.nAtoms() << " no matching atoms "
            << queryAtoms.size();
    REPORT(Reporter::TRACE) << "target molecule no of atoms "
            << targetMolecule.nAtoms() << " no matching atoms "
            << targetAtoms.size();

    // create the adjacency matrices

    size_t querySize = queryAtoms.size();
    queryAtomConnections = std::vector<size_t>(querySize);
    queryAdjacencyMatrix = make_unique<Array2D<bool>>(querySize, querySize);
    buildAdjacencyMatrix(queryMolecule, queryAtoms, *queryAdjacencyMatrix,
            queryAtomConnections);
    REPORT(Reporter::TRACE) << "Query adjacency matrix "
            << *queryAdjacencyMatrix;
    REPORT(Reporter::TRACE) << "Query connections  "
            << collectionToString(queryAtomConnections);

    size_t targetSize = targetAtoms.size();
    targetAtomConnections = std::vector<size_t>(targetSize);
    targetAdjacencyMatrix = make_unique<Array2D<bool>>(targetSize, targetSize);
    buildAdjacencyMatrix(targetMolecule, targetAtoms, *targetAdjacencyMatrix,
            targetAtomConnections);
    REPORT(Reporter::TRACE) << "Target adjacency matrix "
            << *targetAdjacencyMatrix;
    REPORT(Reporter::TRACE) << "Target connections  "
            << collectionToString(targetAtomConnections);

    // Find all isomorphisms
    Ullmann<MolComparer> ullmann(*this);
    auto nIsomorphisms = ullmann.doUllman();

    REPORT(Reporter::DETAIL) << "found " << nIsomorphisms
            << " structure isomorphisms";
    return nIsomorphisms;
}

std::vector<const Atom *> MolComparer::findAtoms(
        const Molecule & molecule) const {
    auto map = [] (const unique_ptr<Atom> & atom) -> const Atom * {
        return atom.get();
    };
    auto atoms = mapToNewList<unique_ptr<Atom>, const Atom *>(
            molecule.getAtoms(), map);
    if (heavyAtomOnly) {
        auto filter = [](const Atom * atom)->bool {
            return AtomType::isHeavy(atom->getAtomType().getType());
        };
        atoms = filterList<const Atom *>(atoms, filter);
    }
    return atoms;
}

bool MolComparer::sameNodeType(const size_t queryNodeNo,
        const size_t targetNodeNo) const {
    if (matchTypes) {
        auto queryAtom = queryAtoms.at(queryNodeNo);
        auto targetAtom = targetAtoms.at(targetNodeNo);
        auto queryType = queryAtom->getAtomType().getType();
        auto targetType = targetAtom->getAtomType().getType();
        if (matchElementalTypes) {
            queryType = AtomType::getElementalType(queryType);
            targetType = AtomType::getElementalType(targetType);
        }
        if (queryType != targetType) {
            return false;
        }
    }
    auto nQueryConnections = queryAtomConnections.at(queryNodeNo);
    auto nTargetConnections = targetAtomConnections.at(targetNodeNo);
    if (subgraph) {
        if (nQueryConnections > nTargetConnections) {
            return false;
        }
    } else {
        if (nQueryConnections != nTargetConnections) {
            return false;
        }
    }
    return true;
}

bool MolComparer::sameEdgeType(const size_t queryNodeNo1,
        const size_t queryNodeNo2, const size_t targetNodeNo1,
        const size_t targetNodeNo2) const {
    return true;
}

// functions required by the Ullman algorithm

size_t MolComparer::getSizeOfTarget() {
    return targetAtoms.size();
}

size_t MolComparer::getSizeOfQuery() {
    return queryAtoms.size();
}

bool MolComparer::isAdjacentInTarget(const size_t targetNodeNo1,
        const size_t targetNodeNo2) {
    return targetAdjacencyMatrix->get(targetNodeNo1, targetNodeNo2);
}

bool MolComparer::isAdjacentInQuery(const size_t queryNodeNo1,
        const size_t queryNodeNo2) {
    return queryAdjacencyMatrix->get(queryNodeNo1, queryNodeNo2);
}

void MolComparer::callback(const std::vector<size_t> & queryIdsToTargetIds) {
    nIsomorphisms++;
    callbackFunction(queryIdsToTargetIds);
}

void MolComparer::defaultCallbackFunction(
        const std::vector<size_t> & queryIdsToTargetIds) {
    cout << "Isomorphism " << nIsomorphisms << endl;
    for (size_t queryId = 0; queryId < queryIdsToTargetIds.size(); queryId++) {
        auto targetId = queryIdsToTargetIds.at(queryId);
        auto queryAtom = queryAtoms.at(queryId);
        auto targetAtom = targetAtoms.at(targetId);
        cout << "Query atom " << queryAtom->info() << " maps to target atom "
                << targetAtom->info() << endl;
    }
}

unique_ptr<MolComparer> MolComparer::createIsomorphismComparer(
        const Molecule &molecule) {
    auto comparer = make_unique<MolComparer>(molecule, molecule);
    comparer->setHeavyAtomOnly(true);
    comparer->setSubgraph(false);
    return comparer;
}

} /* namespace GarethMol */
