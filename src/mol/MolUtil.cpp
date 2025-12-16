//
// Created by gareth on 12/1/22.
//

#include "MolUtil.h"
#include <ForceField/UFF/TorsionAngle.h>

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

    double torsionAngle(const RDGeom::Point3D &p1, const RDGeom::Point3D &p2,
                        const RDGeom::Point3D &p3, const RDGeom::Point3D &p4) {

        RDGeom::Point3D r1 = p1 - p2, r2 = p3 - p2, r3 = p2 - p3, r4 = p4 - p3;
        RDGeom::Point3D t1 = r1.crossProduct(r2);
        RDGeom::Point3D t2 = r3.crossProduct(r4);
        double d1 = t1.length(), d2 = t2.length();
        double cosPhi = t1.dotProduct(t2) / (d1 * d2);
        if (cosPhi > 1.0) {
            cosPhi = 1.0;
        } else if (cosPhi < -1.0){
            cosPhi = -1.0;
        }
        auto tvp = t1.crossProduct(t2);
        auto sp2 = tvp.dotProduct(r3);
        double phi = acos(cosPhi);
        return sp2 < 0 ? phi : -phi;
    }

    double squareDistance(const RDGeom::Point3D &p1, const RDGeom::Point3D &p2) {
        auto diff = p2 - p1;
        return diff.lengthSq();
    }

    int countSingleBonds(const Atom *atom) {
        const auto &mol = atom->getOwningMol();
        int count = 0;
        for (const auto &bond: mol.atomBonds(atom)) {
            auto bondType = bond->getBondType();
            if (bondType == Bond::BondType::SINGLE) {
                count++;
            }
        }
        return count;
    }

    bool likelySp2(const Atom *atom) {
        const auto &mol = atom->getOwningMol();
        for (const auto &bond: mol.atomBonds(atom)) {
            auto bondType = bond->getBondType();
            if (bondType == Bond::BondType::DOUBLE || bondType == Bond::BondType::AROMATIC) {
                return true;
            }
        }
        return false;
    }

    bool hasTripleBonds(const Atom *atom) {
        const auto &mol = atom->getOwningMol();
        for (const auto &bond: mol.atomBonds(atom)) {
            auto bondType = bond->getBondType();
            if (bondType == Bond::BondType::TRIPLE) {
                return true;
            }
        }
        return false;
    }
}
