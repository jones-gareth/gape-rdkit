//
// Created by jones on 11/25/2022.
//

#ifndef GAPE_GAPEAPP_H
#define GAPE_GAPEAPP_H

#include <rapidjson/document.h>
#include <memory>
#include <vector>

namespace Gape {

    struct GapeSettings {
        bool flattenBonds = false;
        bool flipAmideBonds = false;
        bool solvateStructures = false;


        GapeSettings() = default;

        GapeSettings(const GapeSettings &) = delete;

        GapeSettings &operator=(GapeSettings &) = delete;

        GapeSettings(GapeSettings &&) = default;
    };

    struct SolvationRule;
    struct HydrogenBondingType;

    class GapeApp {
    public:

        explicit GapeApp(const std::string &configFile = "");

        GapeApp(const GapeApp &) = delete;

        GapeApp &operator=(GapeApp &) = delete;

        [[nodiscard]] const GapeSettings &getGapeSettings() const { return gapeSettings; }

        [[nodiscard]] const std::vector<std::shared_ptr<const SolvationRule>>
        &getSolvationRules() const { return solvationRules; }

        [[nodiscard]] const std::vector<std::shared_ptr<const HydrogenBondingType>>
        &getHydrogenBondingTypes() const {return hydrogenBondingTypes; }

    private:
        GapeSettings gapeSettings;
        std::vector<std::shared_ptr<const SolvationRule>> solvationRules;
        std::vector<std::shared_ptr<const HydrogenBondingType>> hydrogenBondingTypes;
    };
}

#endif //GAPE_GAPEAPP_H
