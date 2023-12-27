//
// Created by gareth on 10/18/22.
//

#pragma once

#include <memory>
#include <vector>

#include "SuperpositionMolecule.h"

namespace Gape {
    class Superposition {
        const static int MAX_MOLS;
        SuperpositionMolecule *baseMolecule = nullptr, *fittingMolecule = nullptr;

        std::vector<SuperpositionMolPtr> molecules;

        size_t numberMoleculesMapped = 0U;
        std::vector<int> integerEntryPoints;
        std::vector<int> binaryEntryPoints;
        std::vector<int> integerStringRanges;
        int integerStringLength = -1;
        int binaryStringLength = -1;
        int baseMoleculeNumber = -1;
        int fittingMoleculeNumber = -1;

    public:
        Superposition(const std::vector<std::shared_ptr<SuperpositionMolecule>>& m,
                      const GapeSettings& s);

        void setupMolecules();

        void findBaseMolecule();

        [[nodiscard]] int getIntegerStringLength() const { return integerStringLength; }
        [[nodiscard]] int getBinaryStringLength() const { return binaryStringLength; }
        [[nodiscard]] std::vector<int> getIntegerStringRanges() const { return integerStringRanges; }
        [[nodiscard]] const std::vector<std::shared_ptr<SuperpositionMolecule>>& getMolecules() const {return molecules; }

        const GapeSettings& settings;

    private:
        [[nodiscard]] SuperpositionMolecule* leastFlexibleMolecule() const;

        [[nodiscard]] SuperpositionMolecule* mostFlexibleMolecule() const;

        [[nodiscard]] SuperpositionMolecule* mostFeaturedMolecule() const;

        [[nodiscard]] SuperpositionMolecule* mostFeaturedRigidMolecule() const;

        [[nodiscard]] SuperpositionMolecule* mostActiveRigidMolecule() const;

        [[nodiscard]] SuperpositionMolecule* mostActiveMolecule() const;

        void setWeights() const;
    };
} // GapeApp
