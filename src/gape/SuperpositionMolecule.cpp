//
// Created by gareth on 10/18/22.
//

// #include <GraphMol/ForceFieldHelpers/MMFF/MMFF.h>

#include "SuperpositionMolecule.h"
#include "util/Reporter.h"
#include "mol/Solvate.h"
#include "RotatableBond.h"
#include <GraphMol/MolOps.h>
#include <GraphMol/DistGeomHelpers/Embedder.h>
#include <GraphMol/ForceFieldHelpers/MMFF/MMFF.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <util/GaussianList.h>

#include "mol/HydrogenBondingType.h"
#include "mol/AcceptorAtomFeature.h"
#include "mol/AromaticRingFeature.h"
#include "mol/DonorHydrogenFeature.h"
#include "mol/HydrophobicFeature.h"
#include "mol/PartialCharge.h"

using namespace RDKit;

namespace Gape {
    double VdwInfo::vdwEnergy(const RDKit::Conformer &conformer, double cutoffSqr) const {
        const auto &point1 = conformer.getAtomPos(index0);
        const auto &point2 = conformer.getAtomPos(index1);
        const auto diff = point1 - point2;
        double sqrDist = diff.lengthSq();
        if (sqrDist >= cutoffSqr) {
            return .0;
        }

        const auto dist = sqrt(sqrDist);
        const auto energy = ForceFields::MMFF::Utils::calcVdWEnergy(dist, mmffVdw.R_ij_star, mmffVdw.epsilon);
        return energy;
    }

    SuperpositionMolecule::SuperpositionMolecule(const ROMol &inputMol, const GapeSettings &settings) : settings(
        settings) {
        mol = inputMol;
        MolOps::addHs(mol, false, true);
        MolOps::findSSSR(mol);

        mmffMolProperties = new MMFF::MMFFMolProperties(mol);
        assert(mmffMolProperties->isValid());
        if (mol.getNumConformers() == 1) {
            referenceConformer = mol.getConformer();
        } else {
            assert(mol.getNumConformers() == 0);
        }
        auto name = mol.getProp<std::string>("_Name");
        REPORT(Reporter::DEBUG) << "mol " << name << " mum atoms " << mol.getNumAtoms() << " num bonds " << mol.
getNumBonds();
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
        if (settings.getGapeParameters().solvateStructures) {
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
                if (canFlatten && settings.getGapeParameters().flattenBonds) {
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
        return (bond.getBeginAtom()->getTotalDegree() == 1 || bond.getEndAtom()->getTotalDegree() == 1);
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
        bool flipAmideBonds = settings.getGapeParameters().flipAmideBonds;
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

    bool SuperpositionMolecule::isNitroOxygen(const Atom &atom) const {
        if (atom.getAtomicNum() != 8 || atom.getDegree() != 1) {
            return false;
        }

        const auto neighbor = *mol.atomNeighbors(&atom).begin();
        return isNitroNitrogen(*neighbor);
    }

    bool SuperpositionMolecule::isNitroNitrogen(const Atom &atom) const {
        if (atom.getAtomicNum() != 7 || atom.getDegree() != 3) {
            return false;
        }

        int numOxy = 0;
        for (const auto &neighbor: mol.atomNeighbors(&atom)) {
            if (neighbor->getAtomicNum() == 8 && neighbor->getDegree() == 1) {
                numOxy++;
            }
        }

        return numOxy == 2;
    }

    bool SuperpositionMolecule::isCarboxylateOxygen(const Atom &atom) const {
        if (atom.getAtomicNum() != 8 || atom.getDegree() != 1) {
            return false;
        }
        const auto neighbor = *mol.atomNeighbors(&atom).begin();
        return isCarboxylateCarbon(*neighbor);
    }

    bool SuperpositionMolecule::isCarboxylateCarbon(const Atom &atom) const {
        if (atom.getAtomicNum() != 6 || atom.getDegree() != 3) {
            return false;
        }

        int numOxy = 0;
        for (const auto &neighbor: mol.atomNeighbors(&atom)) {
            if (neighbor->getAtomicNum() == 8 && neighbor->getDegree() == 1) {
                numOxy++;
            }
        }

        return numOxy == 2;
    }


    void SuperpositionMolecule::findDonorsAndAcceptors() {
        donors = findHydrogenBondDonors(settings.getHydrogenBondingTypes(), mol);
        acceptors = findHydrogenBondAcceptors(settings.getHydrogenBondingTypes(), mol);
    }

    void SuperpositionMolecule::findCharges() {
        findPartialCharges(settings.getPartialCharges(), mol);
    }

    void SuperpositionMolecule::findFeatures() {
        features[FeatureType::DonorInteractionPoint] = DonorHydrogenFeature::findDonorHydrogens(this);
        features[FeatureType::AcceptorAtomFeature] = AcceptorAtomFeature::findAcceptorAtoms(this);
        features[FeatureType::HydrophobicAtom] = HydrophobicFeature::findHydrophobicFeatures(this);
        features[FeatureType::AromaticRing] = AromaticRingFeature::findAromaticRings(this);

        for (const auto &[featureType, featuresForType]: features) {
            append(allFeatures, featuresForType);
        }

        const std::function filter = [](const shared_ptr<Feature> &f) { return f->isMappingFeature(); };
        allMappingFeatures = filterListToNewList(allFeatures, filter);

        donorHydrogens.clear();
        for (const auto &donorHydrogenFeature: features[FeatureType::DonorInteractionPoint]) {
            donorHydrogens.insert(donorHydrogenFeature->getAtom());
        }

        buildSuperpositionCoordinates();
    }

    std::string SuperpositionMolecule::getName() const {
        return getMol().getProp<std::string>("_Name");
    }

    size_t SuperpositionMolecule::numberFeatures() const {
        return allFeatures.size();
    }

    size_t SuperpositionMolecule::numberMappingFeatures() const {
        return allMappingFeatures.size();
    }

    int SuperpositionMolecule::conformationalBitLen() const {
        int nBits = 0;
        for (const auto &rotatableBond: rotatableBonds) {
            switch (rotatableBond->getRotatableBondType()) {
                case RotatableBondType::Flip:
                    nBits += 1;
                    break;
                case RotatableBondType::Full:
                    nBits += 8;
                    break;
                default:
                    ;
            }
        }
        return nBits;
    }

    void SuperpositionMolecule::buildSuperpositionCoordinates() {
        superpositionCoordinates = std::make_unique<SuperpositionCoordinates>(referenceConformer);
        for (const auto &feature: allFeatures) {
            feature->calculateCoordinates(*superpositionCoordinates);
        }
    }

    double SuperpositionMolecule::calculateConformationalEnergy(const RDKit::Conformer &conformer) const {
        double torsionalEnergy = .0;
        for (const auto &rotatableBond: rotatableBonds) {
            torsionalEnergy += rotatableBond->rotatableBondEnergy(conformer);
        }
        double vdwEnergy = .0;
        const double cutoff = settings.getGapeParameters().vdwCutoff;
        const double cutoffSqr = cutoff * cutoff;
        for (const auto &vdwInfo: pairsToCheck) {
            // skip contacts between donor hydrogens and acceptors
            const auto &atom1 = mol.getAtomWithIdx(vdwInfo.index0);
            const auto &atom2 = mol.getAtomWithIdx(vdwInfo.index1);
            if (atom1->getAtomicNum() == 1 && atom2->getAtomicNum() > 1) {
                if (donorHydrogens.find(atom1) != donorHydrogens.end() &&
                    acceptors.find(atom2) != acceptors.end()) {
                    continue;
                }
            } else if (atom2->getAtomicNum() == 1 && atom1->getAtomicNum() > 1) {
                if (donorHydrogens.find(atom2) != donorHydrogens.end() &&
                    acceptors.find(atom1) != acceptors.end()) {
                    continue;
                }
            }
            vdwEnergy += vdwInfo.vdwEnergy(conformer, cutoffSqr);
        }

        return torsionalEnergy + vdwEnergy;
    }

    /**
     * Returns common molecular volume using atomic Gaussians of hydrophobic features.
     *
     * @param mol
     * @return
     */
    double SuperpositionMolecule::gaussianIntegral(const SuperpositionMolecule &otherMolecule, const Conformer &conformer,
                                                   const Conformer &otherConformer) const {
        const GaussianList gaussians(*this, conformer);
        const GaussianList otherGaussians(otherMolecule, otherConformer);

        auto volume = gaussians.overlapVolume(otherGaussians);
        REPORT(Reporter::DEBUG) << "1st order Gaussian integral: " << volume;

        // Java returned here, higher order contributions were not working
        if (true) return volume;

        GaussianList overlay = gaussians.intersection(otherGaussians);
        volume = overlay.volume();
        REPORT(Reporter::DEBUG) << "1st order Gaussian integral (2): " << volume;
        // Subtract 2 order intersections
        overlay = overlay.intersection();
        double diff = overlay.volume();
        volume -= diff;
        REPORT(Reporter::DEBUG) << "2nd order Gaussian Integral " << diff << " new vol " << volume;

        // Add 3 order intersections
        overlay = overlay.intersection();
        diff = overlay.volume();
        volume += diff;
        REPORT(Reporter::DEBUG) << "3nd order Gaussian Integral " << diff << " new vol " << volume;

        return volume;
    }
} // namespace GAPE
