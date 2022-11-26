//
// Created by gareth on 10/18/22.
//

// #include <GraphMol/ForceFieldHelpers/MMFF/MMFF.h>

#include "SuperpositionMolecule.h"
#include <GraphMol/MolOps.h>
#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <GraphMol/ForceFieldHelpers/MMFF/MMFF.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <ForceField/MMFF/Contribs.h>

using namespace RDKit;

namespace Gape {

    SuperpositionMolecule::SuperpositionMolecule(const ROMol &inputMol, const GapeSettings &settings) : settings(settings) {
        mol = inputMol;
        // TODO solvate
        MolOps::addHs(mol);
        MolOps::findSSSR(mol);
        DGeomHelpers::EmbedParameters embedParameters;
        auto confId = EmbedMolecule(mol, embedParameters);
        MMFF::MMFFMolProperties mmffMolProperties(mol);
        assert(mmffMolProperties.isValid());
        forceField = MMFF::constructForceField(mol, 1000);
        ForceFieldsHelper::OptimizeMolecule(*forceField);

        std::vector<boost::shared_ptr<const MMFF::TorsionAngleContrib>> torsions;
        std::vector<boost::shared_ptr<const MMFF::VdWContrib>> vdwPairs;

        for (auto contrib: forceField->contribs()) {
            const ForceFields::ForceFieldContrib *ptr = contrib.get();
            auto *torsion = dynamic_cast<const MMFF::TorsionAngleContrib *>(ptr);
            if (torsion != nullptr) {
                torsions.push_back(static_cast<boost::shared_ptr<const MMFF::TorsionAngleContrib>>(torsion));
                // ForceFields::MMFF::Utils::calcTorsionEnergy
                continue;
            }
            auto *vdw = dynamic_cast<const MMFF::VdWContrib *>(ptr);
            if (vdw != nullptr) {
                vdwPairs.push_back(static_cast<boost::shared_ptr<const MMFF::VdWContrib>>(vdw));
            }
        }

        for (const auto bond: mol.bonds()) {
            bool canFlatten = false;
            isRotatableBond(*bond, true, canFlatten);
        }
        std::cerr << "mol mum atoms " << mol.getNumAtoms() << " num bonds " << mol.getNumBonds() << " num torsions "
                  << torsions.size() << " num vdw " << vdwPairs.size() << std::endl;
    }

    SuperpositionMolecule::~SuperpositionMolecule() {
        delete forceField;
    }

    std::string SuperpositionMolecule::ToMolBlock() const {
        return MolToMolBlock(mol);
    }

    void SuperpositionMolecule::findFreelyRotatableBonds() {


    }

    bool SuperpositionMolecule::isO2(const Atom &atom) const {
        if (atom.getAtomicNum() != 8 || atom.getDegree() != 1 || atom.getHybridization() != Atom::SP2) {
            return false;
        }
        // Don't think this is needed
        const auto bond = *mol.atomBonds(&atom).begin();
        if (bond->getBondType() == Bond::BondType::DOUBLE) {
            return true;
        }
        return false;
    }

    bool SuperpositionMolecule::isO3(const RDKit::Atom &atom) const {
        if (atom.getAtomicNum() != 8 || atom.getTotalDegree() != 2 || atom.getHybridization() != Atom::SP3) {
            return false;
        }
        // Don't think this is needed
        for (const auto bond: mol.atomBonds(&atom)) {
            if (bond->getBondType() != Bond::BondType::SINGLE) {
                return false;
            }
        }
        return true;
    }

    bool SuperpositionMolecule::isAmideBond(const RDKit::Bond &bond) const {
        if (bond.getBondType() != Bond::BondType::SINGLE) {
            return false;
        }
        const auto atom1 = bond.getBeginAtom();
        const auto atom2 = bond.getEndAtom();
        Atom *carbon = nullptr, *nitrogen = nullptr;
        if (atom1->getAtomicNum() == 6) carbon = atom1;
        if (atom1->getAtomicNum() == 7) nitrogen = atom1;
        if (atom2->getAtomicNum() == 6) carbon = atom2;
        if (atom2->getAtomicNum() == 6) nitrogen = atom2;
        if (carbon == nullptr || nitrogen == nullptr) {
            return false;
        }
        if (carbon->getTotalDegree() != 3 || nitrogen->getTotalDegree() != 3) {
            return false;
        }
        if (carbon->getIsAromatic() || nitrogen->getIsAromatic()) {
            return false;
        }
        for (const auto neighbor: mol.atomNeighbors(carbon)) {
            if (isO2(*neighbor)) {
                return true;
            }
        }
        return false;
    }

    bool SuperpositionMolecule::isTerminalBond(const RDKit::Bond &bond) {
        return (bond.getEndAtom()->getTotalDegree() == 1 || bond.getEndAtom()->getTotalDegree() == 1);
    }

    bool SuperpositionMolecule::isNpl3Atom(const RDKit::Atom &atom) const {
        if (atom.getAtomicNum() != 7 || atom.getTotalDegree() != 3) {
            return false;
        }
        for (const auto neighbor: mol.atomNeighbors(&atom)) {
            if (neighbor->getAtomicNum() != 7) {
                return false;
            }
        }
        return true;
    }

    bool SuperpositionMolecule::isArginineCarbon(const RDKit::Atom &atom) const {
        if (atom.getAtomicNum() != 6 || atom.getHybridization() != Atom::SP2 || atom.getTotalDegree() != 3) {
            return false;
        }
        for (const auto neighbor: mol.atomNeighbors(&atom)) {
            if (neighbor->getAtomicNum() != 7) {
                return false;
            }
        }
        return true;
    }

    bool SuperpositionMolecule::isSp2Carbon(const Atom &atom) {
        return atom.getAtomicNum() == 6 && atom.getHybridization() == Atom::SP2;
    }

    bool SuperpositionMolecule::atomIsInRing(const RDKit::Atom &atom) const {
        return mol.getRingInfo()->numAtomRings(atom.getIdx());
    }

    bool SuperpositionMolecule::isCOOHCarbon(const RDKit::Atom &atom, Atom *&o2Atom, Atom *&o3Atom) const {
        if (atom.getAtomicNum() != 6 || atom.getHybridization() != Atom::SP2 || atomIsInRing(atom) ||
            atom.getDegree() != 3) {
            return false;
        }
        o3Atom = nullptr, o2Atom = nullptr;
        for (const auto neighbor: mol.atomNeighbors(&atom)) {
            if (isO2(*neighbor)) {
                o2Atom = neighbor;
            } else if (isO3(*neighbor)) {
                o3Atom = neighbor;
            }
        }
        if (o2Atom == nullptr || o3Atom == nullptr) {
            return false;
        }
        return o3Atom->getTotalNumHs(true) == 1;
    }

    RotatableBondType
    SuperpositionMolecule::isRotatableBond(const RDKit::Bond &bond, const bool flipAmideBonds, bool &canFlatten) const {
        canFlatten = false;
        if (bond.getBondType() != Bond::BondType::SINGLE) {
            return RotatableBondType::None;
        }
        if (mol.getRingInfo()->numBondRings(bond.getIdx())) {
            return RotatableBondType::None;
        };
        if (isAmideBond(bond)) {
            canFlatten = true;
            return flipAmideBonds ? RotatableBondType::Flip : RotatableBondType::None;
        }
        if (isTerminalBond(bond)) {
            return RotatableBondType::None;
        }

        auto *atom1 = bond.getBeginAtom();
        auto *atom2 = bond.getEndAtom();
        if (atom1->getAtomicNum() < atom2->getAtomicNum()) {
            std::swap(atom1, atom2);
        }

        // bonds between C.2 and N.pl3 don't normally rotate
        // but they can flip by 180 degrees
        // However, N.pl3 linking rings is unlikely to be able to stay planar
        if (isNpl3Atom(*atom1) && !atomIsInRing(*atom1) && isSp2Carbon(*atom2)) {
            canFlatten = true;
            if (atom1->getTotalNumHs(true) == 2) {
                return RotatableBondType::None;
            }
            return RotatableBondType::Flip;
        }

        // =N-N= is not rotatable
        if (atom1->getAtomicNum() == 7 && atom1->getHybridization() == Atom::SP2 && atom2->getAtomicNum() == 7
            && atom2->getHybridization() == Atom::SP2) {
            return RotatableBondType::None;
        }

        // COOH is planar
        if (Atom *o3Atom = nullptr, *o2Atom = nullptr; isCOOHCarbon(*atom2, o2Atom, o3Atom)) {
            if (o3Atom == atom1) {
                canFlatten = true;
            }
            return RotatableBondType::Flip;
        }

        return RotatableBondType::Full;

    }

} // namespace GAPE
