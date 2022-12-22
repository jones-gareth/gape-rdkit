//
// Created by gareth on 12/1/22.
//

#include "MolUtil.h"

namespace Gape {

    void findAtoms(const Atom *const startAtom, const std::vector<const Atom *> &stopAtoms,
                                          std::set<const Atom *> &foundAtoms) {
        std::vector<const Atom *> atomsToCheck = {startAtom};

        while (!atomsToCheck.empty()) {
            std::vector<const Atom *> nextAtomsToCheck;
            for (const auto atomToCheck: atomsToCheck) {
                if (auto found = std::find(stopAtoms.begin(), stopAtoms.end(), atomToCheck); found != stopAtoms.end()) {
                    continue;
                }
                if (foundAtoms.find(atomToCheck) == foundAtoms.end()) {
                    foundAtoms.insert(atomToCheck);
                    const auto &mol2 = atomToCheck->getOwningMol();
                    for (const auto neighbor: mol2.atomNeighbors(atomToCheck)) {
                        nextAtomsToCheck.push_back(neighbor);
                    }
                }
            }
            atomsToCheck = nextAtomsToCheck;
        }
    }

    void
    findAtoms(const Atom *const startAtom, const Atom *const stopAtom,
              std::set<const Atom *> &foundAtoms) {
        std::vector<const Atom *> stopAtoms = {stopAtom};
        findAtoms(startAtom, stopAtoms, foundAtoms);
    }

}