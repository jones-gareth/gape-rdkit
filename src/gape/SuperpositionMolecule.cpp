//
// Created by gareth on 10/18/22.
//

// #include <GraphMol/ForceFieldHelpers/MMFF/MMFF.h>

#include "SuperpositionMolecule.h"
#include "../util/Reporter.h"
#include "../mol/Solvate.h"
#include "RotatableBond.h"
#include <GraphMol/MolOps.h>
#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <GraphMol/ForceFieldHelpers/MMFF/MMFF.h>
#include <GraphMol/FileParsers/FileParsers.h>

#include "mol/HydrogenBondingType.h"
#include "../mol/AcceptorAtomFeature.h"
#include "../mol/DonorHydrogenFeature.h"
#include "../mol/HydrophobicFeature.h"

using namespace RDKit;

namespace Gape {

    SuperpositionMolecule::SuperpositionMolecule(const ROMol &inputMol, const GapeApp &settings) : settings(settings) {
        mol = inputMol;
        MolOps::addHs(mol, false, true);
        MolOps::findSSSR(mol);

        mmffMolProperties = new MMFF::MMFFMolProperties(mol);
        assert(mmffMolProperties->isValid());
        if (mol.getNumConformers() == 1) {
            referenceConformer = mol.getConformer();
        }else {
            assert(mol.getNumConformers() == 0);
        }
        auto name = mol.getProp<std::string>("_Name");
        REPORT(Reporter::DEBUG) << "mol " << name << " mum atoms " << mol.getNumAtoms() << " num bonds " << mol.getNumBonds();
    }

    void SuperpositionMolecule::generate3D() {
        // DGeomHelpers::EmbedParameters embedParameters;
        // DGeomHelpers::EmbedMolecule(mol, embedParameters);
        DGeomHelpers::EmbedMolecule(mol, 10, 1);
        mmffMolProperties = new MMFF::MMFFMolProperties(mol);
        assert(mmffMolProperties->isValid());
        const auto forceField = MMFF::constructForceField(mol, mmffMolProperties, 1000);
        ForceFieldsHelper::OptimizeMolecule(*forceField);
        delete forceField;
        assert(mol.getNumConformers() == 1);
        referenceConformer = mol.getConformer();
        for (const auto &rotatableBond: rotatableBonds) {
            rotatableBond->setTorsionAngles();
        }
    }

    void SuperpositionMolecule::solvate() {
        if (settings.getGapeSettings().solvateStructures) {
            solvateMolecule(settings.getSolvationRules(), mol);
        }
    }

    void SuperpositionMolecule::findFreelyRotatableBonds() {

        rotatableBonds.clear();
        for (const auto bond: mol.bonds()) {
            bool canFlatten = false;
            if (const auto rotatableBondType = isRotatableBond(*bond, canFlatten); rotatableBondType !=
                                                                                   RotatableBondType::None) {
                const auto rotatableBond = std::make_shared<RotatableBond>(rotatableBondType, bond, this);
                if (canFlatten && settings.getGapeSettings().flattenBonds) {
                    // TODO flatten bond
                }
                rotatableBonds.push_back(rotatableBond);
            }

        }
    }

    void SuperpositionMolecule::findPairsToCheck() {
        pairsToCheck.clear();

        auto neighborMatrix = MMFF::Tools::buildNeighborMatrix(mol);
        auto numAtoms = mol.getNumAtoms();

        for (unsigned int i = 0; i < numAtoms; ++i) {
            for (unsigned int j = i + 1; j < numAtoms; ++j) {
                if (MMFF::Tools::getTwoBitCell(neighborMatrix, MMFF::Tools::twoBitCellPos(numAtoms, i, j)) >=
                    MMFF::Tools::RELATION_1_4) {
                    auto atom0 = mol.getAtomWithIdx(i);
                    auto atom1 = mol.getAtomWithIdx(j);
                    bool separatedByRotor = std::find_if(rotatableBonds.begin(), rotatableBonds.end(),
                                                         [atom0, atom1](const auto &rb) {
                                                             return rb->isSeparatedByBond(atom0, atom1);
                                                         }) != rotatableBonds.end();
                    if (separatedByRotor) {
                        if (ForceFields::MMFF::MMFFVdWRijstarEps mmffVdWConstants;
                                mmffMolProperties->getMMFFVdWParams(i, j, mmffVdWConstants)) {
                            pairsToCheck.emplace_back(i, j, mmffVdWConstants);
                        }
                    }
                }
            }
        }
    }

    SuperpositionMolecule::~SuperpositionMolecule() {
        delete mmffMolProperties;
    }

    std::string SuperpositionMolecule::ToMolBlock() const {
        return MolToMolBlock(mol);
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
        // RDKit will set the hybridization of the OH atom in COOH to SP2- so don't check hybridization here
        if (atom.getAtomicNum() != 8 || atom.getTotalDegree() != 2) {
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
        if (atom2->getAtomicNum() == 7) nitrogen = atom2;
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
            if (isSp2Carbon(*neighbor)) {
                return true;
            }
        }
        return false;
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
    SuperpositionMolecule::isRotatableBond(const RDKit::Bond &bond, bool &canFlatten) const {
        bool flipAmideBonds = settings.getGapeSettings().flipAmideBonds;
        canFlatten = false;
        if (bond.getBondType() != Bond::BondType::SINGLE) {
            return RotatableBondType::None;
        }
        auto *atom1 = bond.getBeginAtom();
        auto *atom2 = bond.getEndAtom();
        auto hyb = atom1->getHybridization();
        if (hyb != Atom::HybridizationType::SP3 && hyb != Atom::HybridizationType::SP2) {
            return RotatableBondType::None;
        }
        hyb = atom2->getHybridization();
        if (hyb != Atom::HybridizationType::SP3 && hyb != Atom::HybridizationType::SP2) {
            return RotatableBondType::None;
        }
        if (atom1->getDegree() == 1 || atom2->getDegree() == 1) {
            return RotatableBondType::None;
        }
        if (mol.getRingInfo()->numBondRings(bond.getIdx())) {
            return RotatableBondType::None;
        }
        if (isAmideBond(bond)) {
            canFlatten = true;
            return flipAmideBonds ? RotatableBondType::Flip : RotatableBondType::None;
        }
        if (isTerminalBond(bond)) {
            return RotatableBondType::None;
        }

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

        if (atom1->getAtomicNum() == 8 && atom2->getAtomicNum() == 6) {
            // COOH is planar
            if (Atom *o3Atom = nullptr, *o2Atom = nullptr; isCOOHCarbon(*atom2, o2Atom, o3Atom)) {
                if (o3Atom == atom1) {
                    canFlatten = true;
                }
                return RotatableBondType::Flip;
            }
        }

        return RotatableBondType::Full;
    }

    void SuperpositionMolecule::setConformer(const Conformer &conformer) {
        assert(mol.getNumConformers() == 1);
        const auto id = mol.getConformer().getId();
        mol.removeConformer(id);
        const auto conformerPtr = new Conformer(conformer);
        mol.addConformer(conformerPtr, true);
        assert(mol.getNumConformers() == 1);
    }

    bool SuperpositionMolecule::isNitroOxygen(const Atom& atom) const
    {
	    if (atom.getAtomicNum() != 8 || atom.getDegree() != 1)
	    {
            return false;
	    }

        const auto neighbor = mol.atomNeighbors(&atom).begin().current;
        return isNitroNitrogen(*neighbor);
    }

	bool SuperpositionMolecule::isNitroNitrogen(const Atom& atom) const
	{
	    if (atom.getAtomicNum() != 7 || atom.getDegree() != 3)
	    {
            return false;
	    }

        int numOxy = 0;
        for (const auto &neighbor: mol.atomNeighbors(&atom))
        {
	        if (neighbor->getAtomicNum() == 8 && neighbor->getDegree() == 1)
	        {
                numOxy++;
	        }
        }

        return numOxy == 2;
	}

    bool SuperpositionMolecule::isCarboxylateOxygen(const Atom& atom) const
    {
	    if (atom.getAtomicNum() != 8 || atom.getDegree() != 1)
	    {
            return false;
	    }
        const auto neighbor = mol.atomNeighbors(&atom).begin().current;
        return isCarboxylateCarbon(*neighbor);
    }

    bool SuperpositionMolecule::isCarboxylateCarbon(const Atom& atom) const
    {
	    if (atom.getAtomicNum() != 6 || atom.getDegree() != 3)
	    {
            return false;
	    }

    	int numOxy = 0;
        for (const auto &neighbor: mol.atomNeighbors(&atom))
        {
	        if (neighbor->getAtomicNum() == 8 && neighbor->getDegree() == 1)
	        {
                numOxy++;
	        }
        }

        return numOxy == 2;

    }


    void SuperpositionMolecule::findDonorsAndAcceptors()
    {
        donors = findHydrogenBondDonors(settings.getHydrogenBondingTypes(), mol);
        acceptors = findHydrogenBondAcceptors(settings.getHydrogenBondingTypes(), mol);
    }

    void SuperpositionMolecule::findFeatures()
    {
        features[FeatureType::DonorInteractionPoint] = DonorHydrogenFeature::findDonorHydrogens(this);
        features[FeatureType::AcceptorAtomFeature] = AcceptorAtomFeature::findAcceptorAtoms(this);
        features[FeatureType::HydrophobicAtom] = HydrophobicFeature::findHydrophobicFeatures(this);
    }


} // namespace GAPE
