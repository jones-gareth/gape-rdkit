//
// Created by jones on 1/15/2026.
//

#include "FreeCorner.h"

namespace Gape {

    const double FreeCorner::CORNER_FLAP_TOL = 0.15;

    bool FreeCorner::checkFreeCorder(const Atom * const atomX) {
        const auto &mol = atomX->getOwningMol();
        const auto ringInfo = mol.getRingInfo();
        const auto atomRings = ringInfo->atomMembers(atomX->getIdx());
        if (atomRings.size() != 1) {
            return false;
        }

        const auto ringIndices = ringInfo->atomRings()[atomRings[0]];
        if (ringIndices.size() < 4) {
            return false;
        }

        const Atom* atomC = nullptr;
        const Atom* atomB = nullptr;
        const Bond* bondBX = nullptr;
        const Bond* bondXC = nullptr;

        for (int idx: ringIndices) {
            auto bond = mol.getBondBetweenAtoms(idx, atomX->getIdx());
            if (bond == nullptr) {
                continue;
            }
            if (bond->getBondType() != Bond::BondType::SINGLE) {
                return false;
            }
            if (atomB == nullptr) {
                atomB = mol.getAtomWithIdx(idx);
                bondBX = bond;
            } else if (atomC == nullptr) {
                atomC = mol.getAtomWithIdx(idx);
                bondXC = bond;
            }
        }

        const Atom* atomA = nullptr;
        const Atom* atomD = nullptr;
        const Bond* bondAB = nullptr;
        const Bond* bondCD = nullptr;
        for (int idx: ringIndices) {
            
        }


    }
} // Gape