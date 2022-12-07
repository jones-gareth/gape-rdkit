//
// Created by jones on 11/25/2022.
//

#include <string>
#include "GapeApp.h"
#include "rapidjson/document.h"
#include "rapidjson/writer.h"
#include "rapidjson/stringbuffer.h"
#include "../mol/Solvate.h"

namespace Gape {
    std::string defaultConfigStr = R"JSON(
{
    "molecule": {
        // Flatten bonds when possible
        "flattenBonds": true,
        // Allow amide bonds to flip
        "flipAmideBonds": true,
        // Perform common ionizations at physiological Ph
        "solvateStructures": true
    },
    "solvationRules": [
        {
            "name": "Primary Aliphatic Amine",
            "type": "base",
            "smarts": "[N^3X3H2]-[C^3X4]",
            "pKa": 10.0
        },
        {
            "name": "Secondary Aliphatic Amine",
            "type": "base",
            "smarts": "[N^3X3H](-[C^3X4])-[C^3X4]",
            "pKa": 10.0
        },
        {
            "name": "Tertiary Aliphatic Amine",
            "type": "base",
            "smarts": "[N^3X3](-[C^3X4])(-[C^3X4])-[C^3X4]",
            "pKa": 10.0
        },
        {
            "name": "Guanidinium",
            "type": "base",
            "smarts": "[N^2X2H]=[C^2X3](-N)-N",
            "pKa": 14.4
        },
        {
            "name": "Amidine 1",
            "type": "base",
            "smarts": "[N^2X2H]=[C^2X3]-N",
            "pKa": 12.4
        },
        {
            "name": "Amidine 2",
            "type": "base",
            "smarts": "[N^2X2]=[C^2X3]-[NX3]",
            "pKa": 12.4
        },
        {
            "name": "Ortho Aminopryidine",
            "type": "base",
            "smarts": "[nX2]1:c(-[NH2]):c:[c,n]:c:c:1",
            "pKa": 10.0
        },
        {
            "name": "Para Aminopryidine",
            "type": "base",
            "smarts": "[nX2]1:c:[c,n]:c(-[NH2]):c:c:1",
            "pKa": 10.0
        },
        {
            "name": "Carboxyl",
            "type": "acid",
            "smarts": "[OX2H]-C=O",
            "pKa": 4.8
        },
        {
            "name": "Nitro",
            "type": "acid",
            "smarts": "[OX2H]-N=O",
            "pKa": 11
        },
        {
            "name": "Tetrazole tautomer 1",
            "type": "acid",
            "smarts": "[nHX3]1:n:n:c:n:1",
            "pKa": 3.0
        },
        {
            "name": "Tetrazole tautomer 2",
            "type": "acid",
            "smarts": "[nHX3]1:n:n:n:c:1",
            "pKa": 3.0
        },
        /*
        {"name": "Imidazole 1",
        "type": "base",
        "smarts": "[nX2]1:c:[nX3H]:c:c:1",
        "pKa": 6},
        {"name": "Imidazole 2",
        "type": "acid",
        "smarts": [nX2]1:c:[nX3](-C):c:c:@1",
        "pKa": 6},
        */
        // correct protonated purine from imidazole rule
        {
            "name": "Protonated Purine",
            "type": "acid",
            "smarts": "[nX3H]1:c:n:c2:c:1:n:c:n:c:2",
            "pKa": 1.0
        },
        {
            "name": "Acylsulponamide",
            "type": "acid",
            "smarts": "[NX3H](-C=O)-S(=O)=O",
            "pKa": 4.8
        },
        // pKa unknown, but assumed to be about 7
        {
            "name": "Wyeth heterocycle",
            "type": "acid",
            "smarts": "[#7X3H]1[#6](=O)[#8][#7][#6](=O)1",
            "pKa": 7.0
        },
        {
            "name": "Phosphate",
            "type": "acid",
            "smarts": "[OX2H]-P=O",
            "pKa": 6.0
        },
        {
            "name": "Sulphonic acid",
            "type": "acid",
            "smarts": "[OX2H]-S(=O)=O",
            "pKa": 1.0
        }
    ]
}
    )JSON";

    GapeApp::GapeApp(const std::string &configFile) {
        rapidjson::Document d;
        d.Parse<rapidjson::kParseCommentsFlag>(defaultConfigStr.c_str());

        gapeSettings.flattenBonds = d["molecule"]["flattenBonds"].GetBool();
        gapeSettings.flipAmideBonds = d["molecule"]["flipAmideBonds"].GetBool();
        gapeSettings.solvateStructures = d["molecule"]["solvateStructures"].GetBool();

        auto const &jsonSolvationRules = d["solvationRules"];
        solvationRules.clear();
        solvationRules.reserve(jsonSolvationRules.Size());
        for (auto &jsonRule: jsonSolvationRules.GetArray()) {
            std::string name(jsonRule["name"].GetString());
            std::string typeName(jsonRule["type"].GetString());
            assert(typeName == "acid" || typeName == "base");
            std::string smarts(jsonRule["smarts"].GetString());
            auto pKa = jsonRule["pKa"].GetDouble();
            SolvationType solvationType = typeName == "acid" ? SolvationType::Acid : SolvationType::Base;
            auto rule = std::make_shared<const SolvationRule>(solvationType, name, smarts, pKa);
            solvationRules.push_back(rule);
        }

    }
}