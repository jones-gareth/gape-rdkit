//
// Created by gareth on 10/18/22.
//

#pragma once

#include <GraphMol/GraphMol.h>
#include <ForceField/MMFF/Params.h>
#include <util/ConnectedGraphFinder.h>

#include "GapeSettings.h"
#include "mol/Feature.h"
#include "SuperpositionCoordinates.h"


using namespace RDKit;

namespace ForceFields {
    class ForceField;
}


namespace RDKit::MMFF {
    class MMFFMolProperties;
}

namespace Gape {
    enum class RotatableBondType {
        None,
        Flip,
        Full
    };

    class RotatableBond;

    class VdwInfo {
    public:
        const unsigned int index0, index1;
        const ForceFields::MMFF::MMFFVdWRijstarEps mmffVdw;

        VdwInfo(unsigned int index0, unsigned int index1,
                ForceFields::MMFF::MMFFVdWRijstarEps mmffVdw) : index0(index0), index1(index1), mmffVdw(mmffVdw) {
        }
    };

    class SuperpositionMolecule {
    public:
        explicit SuperpositionMolecule(const ROMol& mol, const GapeSettings& settings);

        virtual ~SuperpositionMolecule();

        SuperpositionMolecule(const SuperpositionMolecule&) = delete;

        SuperpositionMolecule& operator=(SuperpositionMolecule&) = delete;

        std::string ToMolBlock() const;

        const bool& isRigid() const { return rigid; }

        const bool& isFixed() const { return fixed; }

        const double& activity() const { return act; }

        const double& getWeight() const { return weight; }

        const RWMol& getMol() const { return mol; }

        MMFF::MMFFMolProperties* getMMFFMolProperties() const { return mmffMolProperties; }

        void findFreelyRotatableBonds();

        const std::vector<std::shared_ptr<RotatableBond>>& getRotatableBonds() const { return rotatableBonds; }

        int conformationalBitLen() const;

        void findPairsToCheck();

        void generate3D();

        void findFeatures();

        void solvate();

        void setConformer(const Conformer& conformer);

        [[nodiscard]] const Conformer& getReferenceConformer() const { return referenceConformer; }

        [[nodiscard]] const SuperpositionCoordinates& getReferenceCoordinates() const { return *superpositionCoordinates; }

        void findDonorsAndAcceptors();

        void findCharges();

        const std::map<const Atom *, std::shared_ptr<const HydrogenBondingType>>& getDonors() const { return donors; }

        const std::map<const Atom *, std::shared_ptr<const HydrogenBondingType>>& getAcceptors() const {
            return acceptors;
        }

        /**
         * Check to see if this carbon is in a nitro group.
         *
         * @return
         */
        bool isNitroOxygen(const Atom& atom) const;

        /**
         * Checks to see if this oxygen is in a carboxylate group
         *
         * @return
         */
        bool isCarboxylateOxygen(const Atom& atom) const;

        std::string getName() const;
        size_t numberFeatures() const;
        size_t numberMappingFeatures() const;
        double act = 0.0;
        const std::map<const FeatureType, std::vector<std::shared_ptr<Feature>>>& getFeatures() const {
            return features;
        }
        const std::vector<std::shared_ptr<Feature>>& getAllFeatures() const { return allFeatures; }
        const std::vector<std::shared_ptr<Feature>>& getAllMappingFeatures() const { return allMappingFeatures; }

        void copySuperpositionCoordinates(const SuperpositionCoordinates& otherCoordinates) const;

    private:
        RWMol mol;
        MMFF::MMFFMolProperties* mmffMolProperties;
        const GapeSettings& settings;
        std::vector<std::shared_ptr<RotatableBond>> rotatableBonds;
        std::vector<VdwInfo> pairsToCheck;
        std::map<const Atom *, std::shared_ptr<const HydrogenBondingType>> donors;
        std::map<const Atom *, std::shared_ptr<const HydrogenBondingType>> acceptors;
        Conformer referenceConformer;
        std::unique_ptr<SuperpositionCoordinates> superpositionCoordinates;
        std::map<const FeatureType, std::vector<std::shared_ptr<Feature>>> features;
        std::vector<std::shared_ptr<Feature>> allFeatures;
        std::vector<std::shared_ptr<Feature>> allMappingFeatures;
        bool rigid = false;
        bool fixed = false;
        double weight = 1.0;

        bool isO2(const Atom& atom) const;

        bool isO3(const Atom& atom) const;

        bool isAmideBond(const Bond& bond) const;

        bool isNpl3Atom(const Atom& atom) const;

        static bool isTerminalBond(const Bond& bond);

        bool isArginineCarbon(const Atom& atom) const;

        static bool isSp2Carbon(const Atom& atom);

        bool atomIsInRing(const Atom& atom) const;

        bool isCOOHCarbon(const Atom& atom, Atom*& o2Atom, Atom*& o3Atom) const;

        RotatableBondType isRotatableBond(const Bond& bond, bool& canFlatten) const;

        /**
         * Check to see if this carbon is in a nitro group.
         *
         * @return
         */
        bool isNitroNitrogen(const Atom& atom) const;

        /**
         * Checks to see if this carbon is a carboxylate
         *
         * @return
         */
        bool isCarboxylateCarbon(const Atom& atom) const;

        void buildSuperpositionCoordinates();

        friend class Superposition;
    };

    using SuperpositionMolPtr = std::shared_ptr<SuperpositionMolecule>;
} // namespace Gape
