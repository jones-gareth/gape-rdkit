//
// Created by Gareth Jones on 11/25/2022.
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
        bool scaleFitting = true;
        double startFittingRadius = 3.5;
        double finishFittingRadius = 1.5;
        int numberRebuilds = 24;
		bool ignoreVdwAttractive = true;
        bool ignoreTorsion = false;
        int numberRuns = 10;
        bool guessGaParameters = true;
		int numberIslands = 5;
		int populationSize = 100;
		int numberOperations = 60000;
        double selectionPressure = 1.0001;
        bool useNiches = true;
        int nicheSize = 5;
		double nichingOff = 0.6;
		int crossoverWeight = 95;
		int mutationWeight = 95;
		int migrationWeight = 10;
        double donorHydrogenWeight = 1750;
        double acceptorAtomWeight = 1750;
        double aromaticRingWeight = 2000;
        double volumeWeight = 100;
        double conformationalWeight = 10;
        double constraintWeight = 100;
        double vdwCutoff = 100;

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

