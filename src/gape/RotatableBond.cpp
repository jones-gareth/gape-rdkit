//
// Created by jones on 11/25/2022.
//

#include "RotatableBond.h"
#include "SuperpositionMolecule.h"
#include <ForceField/MMFF/Contribs.h>
#include <mol/MolUtil.h>

namespace Gape {
    RotatableBond::RotatableBond(RotatableBondType rotatableBondType, const Bond *const bond,
                                 SuperpositionMolecule *superpositionMolecule) :
            rotatableBondType(rotatableBondType), molecule(superpositionMolecule) {
        std::set<const Atom *> atoms1, atoms2;
        const auto &mol = molecule->getMol();
        findAtoms(bond->getBeginAtom(), bond->getEndAtom(), atoms1);
        findAtoms(bond->getEndAtom(), bond->getBeginAtom(), atoms2);
        if (atoms1.size() <= atoms2.size()) {
            atom1 = bond->getBeginAtom();
            atom2 = bond->getEndAtom();
            atom1List.assign(atoms1.begin(), atoms1.end());
            atom2List.assign(atoms2.begin(), atoms2.end());
        } else {
            atom1 = bond->getEndAtom();
            atom2 = bond->getBeginAtom();
            atom1List.assign(atoms2.begin(), atoms2.end());
            atom2List.assign(atoms1.begin(), atoms1.end());
        }

        const auto mmffMolProperties = molecule->getMMFFMolProperties();

        for (const auto atom0: mol.atomNeighbors(atom1)) {
            if (atom0 != atom2) {
                for (const auto atom3: mol.atomNeighbors(atom2)) {
                    if (atom3 != atom1) {
                        ForceFields::MMFF::MMFFTor mmffTorsion;
                        const auto index0 = atom0->getIdx();
                        const auto index1 = atom1->getIdx();
                        const auto index2 = atom2->getIdx();
                        const auto index3 = atom3->getIdx();
                        unsigned int torsionType;
                        if (mmffMolProperties->getMMFFTorsionParams(mol, index0, index1, index2, index3, torsionType,
                                                                    mmffTorsion)) {
                            TorsionInfo torsionInfo(index0, index1, index2, index3, mmffTorsion);
                            torsions.push_back(torsionInfo);
                        }
                    }
                }
            }
        }
    }

    RotatableBond::~RotatableBond() {

    }

    bool RotatableBond::isSeparatedByBond(const Atom * const a1, const Atom * const a2) const {
        if (std::find(atom1List.begin(), atom1List.end(), a1) != atom1List.end() &&
            std::find(atom2List.begin(), atom2List.end(), a2) != atom2List.end()) {
            return true;
        }
        if (std::find(atom1List.begin(), atom1List.end(), a2) != atom1List.end() &&
            std::find(atom2List.begin(), atom2List.end(), a1) != atom2List.end()) {
            return true;
        }
        return false;
    }
}