//
// Created by Gareth Jones on 11/25/2022.
//

#include "RotatableBond.h"
#include "SuperpositionMolecule.h"
#include <ForceField/MMFF/TorsionAngle.h>
#include <Geometry/Transform3D.h>
#include "mol/MolUtil.h"
#include "util/TransformOps.h"

namespace Gape {
    double TorsionInfo::torsionEnergy(const RDKit::Conformer &conformer) const {
        const auto &iPoint = conformer.getAtomPos(index0);
        const auto &jPoint = conformer.getAtomPos(index1);
        const auto &kPoint = conformer.getAtomPos(index2);
        const auto &lPoint = conformer.getAtomPos(index3);

        const auto d_V1 = mmffTorsion.V1;
        const auto d_V2 = mmffTorsion.V2;
        const auto d_V3 = mmffTorsion.V3;
        const auto cosPhi =
                ForceFields::MMFF::Utils::calcTorsionCosPhi(iPoint, jPoint, kPoint, lPoint);

        return ForceFields::MMFF::Utils::calcTorsionEnergy(d_V1, d_V2, d_V3, cosPhi);
    }

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
        setTorsionAngles();
    }

    void RotatableBond::setTorsionAngles() {
        auto const &conformer = molecule->getReferenceConformer();
        if (conformer.getNumAtoms() == 0) {
            return;
        }
        for (auto &torsion: torsions) {
            const auto angle = torsionAngle(conformer.getAtomPos(torsion.index0), conformer.getAtomPos(torsion.index1),
                                            conformer.getAtomPos(torsion.index2), conformer.getAtomPos(torsion.index3));
            torsion.referenceAngle = angle;
        }
    }

    RotatableBond::~RotatableBond() {
    }

    bool RotatableBond::isSeparatedByBond(const Atom *const a1, const Atom *const a2) const {
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

    void RotatableBond::rotateBond(double angle, Conformer &conf) const {
        RDGeom::Transform3D rot;
        assert(conf.getNumAtoms() == molecule->getMol().getNumAtoms());
        determineRotation(conf.getAtomPos(atom1->getIdx()), conf.getAtomPos(atom2->getIdx()), angle, rot);
        for (const auto atom: atom1List) {
            rot.TransformPoint(conf.getAtomPos(atom->getIdx()));
        }
    }

    void RotatableBond::rotateBond(double angle, SuperpositionCoordinates &superpositionCoordinates) const {
        RDGeom::Transform3D rot;
        auto &conf = superpositionCoordinates.getConformer();
        assert(conf.getNumAtoms() == molecule->getMol().getNumAtoms());
        determineRotation(conf.getAtomPos(atom1->getIdx()), conf.getAtomPos(atom2->getIdx()), angle, rot);
        for (const auto atom: atom1List) {
            rot.TransformPoint(conf.getAtomPos(atom->getIdx()));
            for (auto &features: superpositionCoordinates.getFeatureCoordinates()) {
                auto &featureMap = features.second;
                if (auto featureCoordinate = featureMap.find(atom); featureCoordinate != featureMap.end()) {
                    for (auto &coordinate: featureCoordinate->second) {
                        rot.TransformPoint(coordinate);
                    }
                }
            }
        }
    }

    double RotatableBond::rotatableBondEnergy(const RDKit::Conformer &conformer) const {
        double energy = 0;
        for (const auto torsion: torsions) {
            energy += torsion.torsionEnergy(conformer);
        }
        return energy;
    }
}
