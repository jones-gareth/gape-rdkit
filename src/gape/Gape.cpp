//
// Created by jones on 11/25/2022.
//

#include <string>
#include "Gape.h"
#include "rapidjson/document.h"
#include "rapidjson/writer.h"
#include "rapidjson/stringbuffer.h"

namespace Gape {
    std::string defaultConfigStr = R"JSON(
{
   // Flatten bonds when possible
   "flattenBonds": true,
   // Allow amide bonds to flip
   "flipAmideBonds": true,
   // Perform common ionizations at physiological Ph
   "solvateStructures": true
}
    )JSON";

    GapeSettings Gape::readSettings() {
        rapidjson::Document d;
        d.Parse<rapidjson::kParseCommentsFlag>(defaultConfigStr.c_str());

        GapeSettings settings;
        settings.flattenBonds = d["flattenBonds"].GetBool();
        settings.flipAmideBonds = d["flipAmideBonds"].GetBool();
        settings.solvateStructures = d["solvateStructures"].GetBool();

        return settings;
    }
}