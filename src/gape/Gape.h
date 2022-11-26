//
// Created by jones on 11/25/2022.
//

#ifndef GAPE_GAPE_H
#define GAPE_GAPE_H

namespace Gape {

    struct GapeSettings {
        bool flattenBonds;
        bool flipAmideBonds;
        bool solvateStructures;

        GapeSettings() {}
        GapeSettings(const GapeSettings &) = delete;

        GapeSettings &operator=(GapeSettings &) = delete;

        GapeSettings(GapeSettings&&)=default;
    };

    class Gape {
    public:
        const GapeSettings gapeSettings;

        Gape() : gapeSettings(readSettings()) {
        }

        Gape(const Gape &) = delete;

        Gape &operator=(Gape &) = delete;

    private:
        static GapeSettings readSettings();
    };
}

#endif //GAPE_GAPE_H
