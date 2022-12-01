//
// Created by jones on 11/25/2022.
//

#ifndef GAPE_GAPE_H
#define GAPE_GAPE_H

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

    class SolvationRule;

    class Gape {
    public:

        explicit Gape(const std::string &configFile = "");

        Gape(const Gape &) = delete;

        Gape &operator=(Gape &) = delete;

        const GapeSettings &getGapeSettings() const { return gapeSettings; }

        const std::vector<std::shared_ptr<const SolvationRule>> getSolvationRules() const { return solvationRules; }

    private:
        GapeSettings gapeSettings;
        std::vector<std::shared_ptr<const SolvationRule>> solvationRules;
    };
}

#endif //GAPE_GAPE_H
