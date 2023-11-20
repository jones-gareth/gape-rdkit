//
// Created by jones on 11/25/2022.
//

#pragma once

#include <rapidjson/document.h>
#include <memory>
#include <vector>

namespace Gape {

    enum class BaseMoleculeSelection
    {
	    minRotatableBonds,
        maxFeatures,
        maxActivity
    };

    enum class FittingMoleculeSelection
    {
	    minRotatableBonds,
        maxFeatures,
        baseMolecule
    };

    struct GapeParameters {
        bool flattenBonds = false;
        bool flipAmideBonds = false;
        bool solvateStructures = false;
        BaseMoleculeSelection baseMoleculeSelection = BaseMoleculeSelection::minRotatableBonds;
        FittingMoleculeSelection fittingMoleculeSelection = FittingMoleculeSelection::baseMolecule;
        bool useActivities = false;

        GapeParameters() = default;

        GapeParameters(const GapeParameters &) = delete;

        GapeParameters &operator=(GapeParameters &) = delete;

        GapeParameters(GapeParameters &&) = default;
    };

    struct SolvationRule;
    struct HydrogenBondingType;
    struct PartialCharge;

    class GapeSettings {
    public:
        explicit GapeSettings(const std::string &configFile = "");

        GapeSettings(const GapeSettings &) = delete;

        GapeSettings &operator=(GapeSettings &) = delete;

        [[nodiscard]] const GapeParameters &getGapeParameters() const { return gapeParameters; }

        [[nodiscard]] const std::vector<std::shared_ptr<const SolvationRule>>
        &getSolvationRules() const { return solvationRules; }

        [[nodiscard]] const std::vector<std::shared_ptr<const HydrogenBondingType>>
        &getHydrogenBondingTypes() const {return hydrogenBondingTypes; }

        [[nodiscard]] const std::vector<std::shared_ptr<const PartialCharge>>
        &getPartialCharges() const {return partialCharges; }
        
    private:
        GapeParameters gapeParameters;
        std::vector<std::shared_ptr<const SolvationRule>> solvationRules;
        std::vector<std::shared_ptr<const HydrogenBondingType>> hydrogenBondingTypes;
        std::vector<std::shared_ptr<const PartialCharge>> partialCharges;
    };
}

