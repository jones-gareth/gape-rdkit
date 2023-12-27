//
// Created by jones on 11/25/2022.
//

#include <string>
#include "GapeSettings.h"
#include "rapidjson/document.h"
#include "rapidjson/writer.h"
#include "mol/Solvate.h"
#include "mol/HydrogenBondingType.h"
#include "mol/PartialCharge.h"

namespace Gape {
    std::string defaultConfigStr = R"JSON(
{
    "gapeParameters": {
        // Flatten bonds when possible
        "flattenBonds": true,
        // Allow amide bonds to flip
        "flipAmideBonds": true,
        // Perform common ionizations at physiological Ph
        "solvateStructures": true,
		"baseMoleculeSelection": "minRotatableBonds",
		"fittingMoleculeSelection": "baseMolecule",
		"useActivities": false,
		"scaleFitting": true,
		"startFittingRadius": 3.5,
		"finishFittingRadius": 1.5,
		"numberRebuilds": 24,
        "ignoreVdwAttractive": true,
        "ignoreTorsion": false,
		"numberRuns": 10,
        "guessGaParameters": true,
		"numberIslands": 5,
		"populationSize": 100,
		"numberIterations": 60000,
		"selectionPressure": 1.001,
		"useNiches": true,
		"useNiches": true,
        "nicheSize": 5,
		"nichingOff": 0.6,
		"crossoverWeight": 95,
		"mutationWeight": 95,
		"migrationWeight": 10,

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
        "smarts": [nX2]1:c:[nX3](-C):c:c:1",
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
    ],
    "hydrogenBondingTypes": [
        /*  Donor and acceptors from Mills & Dean JCAMD, 10 (1996) 607-622.

            Tabs separate fields

            In the case of cis/trans preferences take the probability mean- if
            only one is present then divide by two

            The sln patterns only support these attributes: f, !h, is= and
            not=.  Z and Y act as Hev but must match the same elemental types,
            additionally X != Y.

            The patterns for amide groups will match N.pl3 atoms.  In GAPE
            these are modeled as planar, not pyramidal.	A correction is
            applied within GAPE.

            If not set weight defaults to number of atoms
        */
        {
            "type": "acceptor",
            "name": "Terminal Phosphate",
            "probability": 0.90,
            "geometry": "cone",
            "smarts": "[OX1]=,-P(=,-[OX1])(-*)=,-[OX1]"
        },
        {
            "type": "acceptor",
            "name": "Phosphinyl",
            "probability": 0.64,
            "geometry": "cone",
            "smarts": "[OX1]=,-P(=,-*)(-*)=,-[OX1]"
        },
        {
            "type": "acceptor",
            "name": "Carboxylate",
            "probability": 0.56,
            "geometry": "dir",
            "smarts": "[OX1]=,-C(-*)=,-[OX1]"
        },
        {
            "type": "acceptor",
            "name": "Asymmetric Het6 N",
            "probability": 0.54,
            "geometry": "dir",
            "smarts": "[#7X2]1~Z~*~*~*~Y~1"
        },
        {
            "type": "acceptor",
            "name": "Acid O",
            "probability": 0.40,
            "geometry": "dir",
            "smarts": "[OX1]=[P,C,S]-[OH]"
        },
        {
            "type": "acceptor",
            "name": "Symmetric Het5 N",
            "probability": 0.48,
            "geometry": "dir",
            "smarts": "[#7X2]1~Z~*~*~Z~1"
        },
        {
            "type": "acceptor",
            "name": "Terminal Sulfate",
            "probability": 0.35,
            "geometry": "cone",
            "smarts": "[OX1]=,-S(=,-[OX1])=,-[OX1]"
        },
        {
            "type": "acceptor",
            "name": "Sulfoxide",
            "probability": 0.40,
            "geometry": "cone",
            "smarts": "[OX1]=S(-*)-*"
        },
        {
            "type": "acceptor",
            "name": "Primary Amine",
            "probability": 0.40,
            "geometry": "dir",
            "smarts": "[NX3H2]-*"
        },
        {
            "type": "acceptor",
            "name": "Asymmetric Het5 N",
            "probability": 0.39,
            "geometry": "dir",
            "smarts": "[#7X2]1~Z~*~*~Y~1"
        },
        {
            "type": "acceptor",
            "name": "Thiocarbonyl",
            "probability": 0.39,
            "geometry": "plane",
            "smarts": "[SX1]=C"
        },
        {
            "type": "acceptor",
            "name": "Hydroxyl",
            "probability": 0.38,
            "geometry": "plane",
            "smarts": "[OX2H]-*"
        },
        {
            "type": "acceptor",
            "name": "Symmetric Het6 N",
            "probability": 0.38,
            "geometry": "dir",
            "smarts": "[#7X2]1~Z~*:~*~*:~Z~1"
        },
        {
            "type": "acceptor",
            "name": "Tertiary Amine",
            "probability": 0.32,
            "geometry": "dir",
            "smarts": "[NX3](*)(*)*"
        },
        {
            "type": "acceptor",
            "name": "Ester O",
            "probability": 0.17,
            "geometry": "plane",
            "smarts": "[OX1]=C(-*)-O-*"
        },
        {
            "type": "acceptor",
            "name": "Carbamate O",
            "probability": 0.20,
            "geometry": "plane",
            "smarts": "[OX1]=C(-O-*)N-*"
        },
        {
            "type": "acceptor",
            "name": "Nitrile",
            "probability": 0.21,
            "geometry": "dir",
            "smarts": "[NX1]#C"
        },
        {
            "type": "acceptor",
            "name": "Imine",
            "probability": 0.19,
            "geometry": "dir",
            "smarts": "[NX2]=C"
        },
        /*
        Dean and Mills have the probability for =O to be 0.18.	 However,  this
        doesn't square with the importance of this group in QSAR studies.
        So, based on empirical observation, I've adjusted the probability to
        that of Thicarbonyl 0.4.
        {"type": "acceptor", "name": "Ketone", "probability": 0.18, "geometry": "plane", "smarts": "[OX1]=C(-*)")},
        */
        {
            "type": "acceptor",
            "name": "Ketone",
            "probability": 0.40,
            "geometry": "plane",
            "smarts": "[OX1]=C(-*)"
        },
        {
            "type": "acceptor",
            "name": "Secondary Amine",
            "probability": 0.18,
            "geometry": "dir",
            "smarts": "[NX3H](-*)-*"
        },
        {
            "type": "acceptor",
            "name": "Phenol OH",
            "probability": 0.13,
            "geometry": "plane",
            "smarts": "[OX2H1]-c1ccccc1"
        },
        {
            "type": "acceptor",
            "name": "Ether",
            "probability": 0.11,
            "geometry": "plane",
            "smarts": "[OX2](-*)-*"
        },
        /*
        In GAPE phenyl NH2 is planar
        {"type": "acceptor", "name": "Phenyl NH2", "probability": 0.08, "geometry": "dir", "smarts": "[NX3H2]-c1ccccc1")},
        */
        // two acceptor geometries for Het5 O, depending on if O is aromatic
        {
            "type": "acceptor",
            "name": "Het5 O",
            "probability": 0.04,
            "geometry": "plane",
            "smarts": "[OX2]1-A-,=A-,=A-,=A-1"
        },
        {
            "type": "acceptor",
            "name": "Het5 O Aromatic",
            "probability": 0.04,
            "geometry": "dir",
            "smarts": "[oX2]1~a~a~a~a~1"
        },
        {
            "type": "acceptor",
            "name": "Nitro O",
            "probability": 0.04,
            "geometry": "dir",
            "smarts": "[OX1]=,-N=,-[OX1]"
        },
        {
            "type": "acceptor",
            "name": "Sulphone",
            "probability": 0.02,
            "geometry": "dir",
            "smarts": "[OX1]=S(=O)-*"
        },
        {
            "type": "donor",
            "name": "Het5 NH",
            "probability": 0.89,
            "smarts": "[#7H]1~*~*~*~*~1"
        },
        {
            "type": "donor",
            "name": "sp2 N+H (1)",
            "probability": 0.82,
            "smarts": "[#7H+](-*)~[#6]~[#7]"
        },
        {
            "type": "donor",
            "name": "sp2 N+H (2)",
            "probability": 0.82,
            "smarts": "[#7H](-*)~[#6]~[#7+]"
        },
        // added charge to these guys, so pyrrole doesn't match
        {
            "type": "donor",
            "name": "sp2 N+H (3)",
            "probability": 0.82,
            "smarts": "[#7H+](:,=*):,-*"
        },
        // Currently, don't support partial charge- but may need this pattern in the future
        /*
        {
            "type": "donor",
            "name": "sp2 N+H (4)",
            "probability": 0.82,
            "smarts": "N[f;fcharge=0.5]H(=:Hev)-:Hev"
        }
        */
        {
            "type": "donor",
            "name": "Acid OH",
            "probability": 0.82,
            "smarts": "[OH]-[C,P,S]=,:O"
        },
        {
            "type": "donor",
            "name": "sp2 N+H2 (1)",
            "probability": 0.75,
            "smarts": "[NH2+]=,:*"
        },
        // sp2 pattern for charged ortho amino pyridine
        {
            "type": "donor",
            "name": "sp2 N+H2 (2)",
            "probability": 0.75,
            "smarts": "[NH2]c1c[c,n]cc[nH+]1"
        },
        // sp2 pattern for charged para amino pyridine
        {
            "type": "donor",
            "name": "sp2 N+H2 (3)",
            "probability": 0.75,
            "smarts": "[NH2]c1[c,n]c[nH+]cc1"
        },
        {
            "type": "donor",
            "name": "sp2 N+H2 (4)",
            "probability": 0.75,
            "smarts": "[#7H2+]~[#6]~[#7]"
        },
        {
            "type": "donor",
            "name": "sp2 N+H2 (5)",
            "probability": 0.75,
            "smarts": "[#7H2]~[#6]~[#7+]"
        },
        {
            "type": "donor",
            "name": "sp3 N+H3",
            "probability": 0.74,
            "smarts": "[NH3+]-*"
        },
        {
            "type": "donor",
            "name": "sp3 N+H2",
            "probability": 0.73,
            "smarts": "[NH2+](-*)-*"
        },
        /*
        Dean and Mills have the probability for N+H to be 0.34.  However this
        doesn't square with the importance of this group in QSAR studies.
        So, based on empirical observation, I've adjusted the probability to
        that of N+H2.
        {"type": "donor", "name": "sp3 N+H", "probability": 0.34, "smarts": "[N+H](*)(*)*"},
        */
        {
            "type": "donor",
            "name": "sp3 N+H",
            "probability": 0.73,
            "smarts": "[N+H](-*)(-*)-*"
        },
        {
            "type": "donor",
            "name": "OH",
            "probability": 0.68,
            "smarts": "[OH]-*"
        },
        {
            "type": "donor",
            "name": "Primary Amide NH2",
            "probability": 0.62,
            "smarts": "[NX3H2]-C=O",
            "weight": 5
        },
        {
            "type": "donor",
            "name": "Secondary Amide NH",
            "probability": 0.54,
            "smarts": "[NX3H](-*)-C=O",
            "weight": 5
        },
        {
            "type": "donor",
            "name": "Phenyl NH",
            "probability": 0.58,
            "smarts": "[NHX3](-*)-c1ccccc1"
        },
        {
            "type": "donor",
            "name": "Phenyl NH2",
            "probability": 0.46,
            "smarts": "[NH2X3]-c1ccccc1"
        },
        {
            "type": "donor",
            "name": "Imine NH",
            "probability": 0.33,
            "smarts": "[NX2H]=*"
        },
        {
            "type": "donor",
            "name": "Phenyl OH",
            "probability": 0.30,
            "smarts": "[OX2H1]-c1ccccc1"
        },
        {
            "type": "donor",
            "name": "Primary Amine NH2",
            "probability": 0.26,
            "smarts": "[NX3H2]-*"
        },
        {
            "type": "donor",
            "name": "Secondary Amine NH",
            "probability": 0.24,
            "smarts": "[NX3H1](-*)-*"
        }
    ], )JSON"
            R"JSON(
	"partialChargeGroups": [
		{
			"name": "Guanidinium 1",
			"formal": 1,
			"partial": 0.333,
			"smarts": "[N^2H2X3]=C(-N)-N",
			"multiple": false
		},
		{
			"name": "Guanidinium 2",
			"formal": 0,
			"partial": 0.333,
			"smarts": "N-C(=[N^2H2X3])-N",
			"multiple": false
		},
		{
			"name": "Amidine 1",
			"formal": 0,
			"partial": 0.5,
			"smarts": "[NX3]-C=[NX3]",
			"multiple": false
		},
		{
			"name": "Amidine 2",
			"formal": 1,
			"partial": 0.5,
			"smarts": "[NX3]=C-[NX3]",
			"multiple": false
		},
		{
			"name": "Ortho Aminopryidine 1",
			"formal": 1,
			"partial": 0.5,
			"smarts": "[nX3H1]1:c(-[NH2]):c:[c,n]:c:c:1",
			"multiple": false
		},
		{
			"name": "Ortho Aminopryidine 2",
			"formal": 0,
			"partial": 0.5,
			"smarts": "[NH2]-c1:[nX3H1]:c:c:[c,n]:c:1",
			"multiple": false
		},
		{
			"name": "Para Aminopryidine 1",
			"formal": 1,
			"partial": 0.5,
			"smarts": "[nX3H1]1:c:[c,n]:c(-[NH2]):c:c:1",
			"multiple": false
		},
		{
			"name": "Para Aminopryidine 2",
			"formal": 0,
			"partial": 0.5,
			"smarts": "[NH2]-c1:c:c:[nX3H1]:c:[c,n]:1",
			"multiple": false
		},
		{
			"name": "Carboxylate 1",
			"formal": -1,
			"partial": -0.5,
			"smarts": "[OX1]-C=O",
			"multiple": false
		},
		{
			"name": "Carboxylate 2",
			"formal": 0,
			"partial": -0.5,
			"smarts": "O=C-[OX1]",
			"multiple": false
		},
		{
			"name": "Hydroxamate 1",
			"formal": 1,
			"partial": -0.5,
			"smarts": "[OX1]-N=O",
			"multiple": false
		},
		{
			"name": "Hydroxamate 2",
			"formal": 0,
			"partial": -0.5,
			"smarts": "O=N-[OX1]",
			"multiple": false
		},
		{
			"name": "Tetrazole 1",
			"formal": 0,
			"partial": -0.25,
			"smarts": "[nX2]1:n:n:c:n:1",
			"multiple": false
		},
		{
			"name": "Tetrazole 2",
			"formal": -1,
			"partial": -0.25,
			"smarts": "n[X2]1:n:n:n:c:1",
			"multiple": false
		},
		{
			"name": "Acylsulponamide",
			"formal": -1,
			"partial": -1.0,
			"smarts": "[NX2H](-C=O)-S(=O)=O",
			"multiple": false
		},
		{
			"name": "Wyeth heterocycle",
			"formal": -1,
			"partial": -1.0,
			"smarts": "[NX2]1-C(=O)-O-N-C(=O)-1",
			"multiple": false
		},
		{
			"name": "Phosphinyl 1",
			"formal": -1,
			"partial": -0.6667,
			"smarts": "[OX1]-P(=O)-[OX1]",
			"multiple": true
		},
		{
			"name": "Phosphinyl 2",
			"formal": 0,
			"partial": -0.6667,
			"smarts": "O=P(-[OX1])-[OX1]",
			"multiple": true
		},
		{
			"name": "Phosphinyl 3",
			"formal": -1,
			"partial": -0.5,
			"smarts": "[OX1]-P=O",
			"multiple": false
		},
		{
			"name": "Phosphinyl 4",
			"formal": 0,
			"partial": -0.5,
			"smarts": "O=P-[OX1]",
			"multiple": false
		},
		{
			"name": "Sulphonyl 1",
			"formal": -1,
			"partial": -0.5,
			"smarts": "[OX1]-S(=O)=O",
			"multiple": false
		},
		{
			"name": "Sulphonyl 2",
			"formal": 0,
			"partial": -0.5,
			"smarts": "O=S(=O)-[OX1]",
			"multiple": false
		},
		{
			"name": "Imidazole 1",
			"formal": 1,
			"partial": 0.5,
			"smarts": "[nX3H1]1:c:[nX3H1]:c:c:1",
			"multiple": false
		},
		{
			"name": "Imidazole 2",
			"formal": 1,
			"partial": 1.0,
			"smarts": "[nX3H1]1:c:[nX3H0]:c:c:1",
			"multiple": false
		},
		{
			"name": "Pyridine",
			"formal": 1,
			"partial": 1.0,
			"smarts": "[nX3H1]1:c:c:c:c:c:1",
			"multiple": false
		},
	    {
			"name": "Charged Amine",
			"formal": 1,
			"partial": 1.0,
			"smarts": "[N^3X4]",
			"multiple": false
		},
	    {
			"name": "N sp2+",
			"formal": 1,
			"partial": 1.0,
			"smarts": "[N^2X3]=*",
			"multiple": false
		}
	]
}
    )JSON";

    GapeSettings::GapeSettings(const std::string& configFile) {
        rapidjson::Document d;
        d.Parse<rapidjson::kParseCommentsFlag>(defaultConfigStr.c_str());

        gapeParameters.flattenBonds = d["gapeParameters"]["flattenBonds"].GetBool();
        gapeParameters.flipAmideBonds = d["gapeParameters"]["flipAmideBonds"].GetBool();
        gapeParameters.solvateStructures = d["gapeParameters"]["solvateStructures"].GetBool();
        std::string baseMolSelection = d["gapeParameters"]["baseMoleculeSelection"].GetString();
        if (baseMolSelection == "minRotatableBonds") {
            gapeParameters.baseMoleculeSelection = BaseMoleculeSelection::minRotatableBonds;
        } else if (baseMolSelection == "maxFeatures") {
            gapeParameters.baseMoleculeSelection = BaseMoleculeSelection::maxFeatures;
        } else {
            REPORT(Reporter::FATAL) << "Unknown base molecule selection " << baseMolSelection;
            throw std::runtime_error("Unknown base molecule selection");
        }
        std::string fittingMolSelection = d["gapeParameters"]["fittingMoleculeSelection"].GetString();
        if (fittingMolSelection == "minRotatableBonds") {
            gapeParameters.fittingMoleculeSelection = FittingMoleculeSelection::minRotatableBonds;
        } else if (fittingMolSelection == "maxFeatures") {
            gapeParameters.fittingMoleculeSelection = FittingMoleculeSelection::maxFeatures;
        } else if (fittingMolSelection == "baseMolecule") {
            gapeParameters.fittingMoleculeSelection = FittingMoleculeSelection::baseMolecule;
        } else {
            REPORT(Reporter::FATAL) << "Unknown fitting molecule selection " << fittingMolSelection;
            throw std::runtime_error("Unknown fitting molecule selection");
        }
        gapeParameters.useActivities = d["gapeParameters"]["useActivities"].GetBool();
        gapeParameters.scaleFitting = d["gapeParameters"]["scaleFitting"].GetBool();
        gapeParameters.startFittingRadius = d["gapeParameters"]["startFittingRadius"].GetDouble();
        gapeParameters.finishFittingRadius = d["gapeParameters"]["finishFittingRadius"].GetDouble();
        gapeParameters.numberRebuilds = d["gapeParameters"]["numberRebuilds"].GetInt();
        gapeParameters.ignoreVdwAttractive = d["gapeParameters"]["ignoreVdwAttractive"].GetBool();
        gapeParameters.ignoreTorsion = d["gapeParameters"]["ignoreTorsion"].GetBool();
        gapeParameters.numberRuns = d["gapeParameters"]["numberRuns"].GetInt();
        gapeParameters.guessGaParameters =d["gapeParameters"]["guessGaParameters"].GetBool();
        gapeParameters.numberIslands =d["gapeParameters"]["numberIslands"].GetInt();
        gapeParameters.populationSize = d["gapeParameters"]["populationsSize"].GetInt();
        gapeParameters.numberIterations = d["gapeParameters"]["numberIterations"].GetInt();
        gapeParameters.selectionPressure = d["gapeParameters"]["selectionPressure"].GetInt();
        gapeParameters.useNiches = d["gapeParameters"]["useNiches"].GetBool();
        gapeParameters.nicheSize = d["gapeParameters"]["nicheSize"].GetInt();
        gapeParameters.nichingOff = d["gapeParameters"]["nichingOff"].GetDouble();
        gapeParameters.crossoverWeight = d["gapeParameters"]["crossoverWeight"].GetInt();
        gapeParameters.mutationWeight = d["gapeParameters"]["mutationWeight"].GetInt();
        gapeParameters.migrationWeight = d["gapeParameters"]["migrationWeight"].GetInt();

        auto const& jsonSolvationRules = d["solvationRules"];
        solvationRules.clear();
        solvationRules.reserve(jsonSolvationRules.Size());
        for (auto& jsonRule: jsonSolvationRules.GetArray()) {
            std::string name(jsonRule["name"].GetString());
            std::string typeName(jsonRule["type"].GetString());
            assert(typeName == "acid" || typeName == "base");
            std::string smarts(jsonRule["smarts"].GetString());
            auto pKa = jsonRule["pKa"].GetDouble();
            SolvationType solvationType = typeName == "acid" ? SolvationType::Acid : SolvationType::Base;
            auto rule = std::make_shared<const SolvationRule>(solvationType, name, smarts, pKa);
            solvationRules.push_back(rule);
        }

        auto const& jsonHydrogenBondingTypes = d["hydrogenBondingTypes"];
        hydrogenBondingTypes.clear();
        hydrogenBondingTypes.reserve(jsonHydrogenBondingTypes.Size());
        for (auto& jsonType: jsonHydrogenBondingTypes.GetArray()) {
            std::string name(jsonType["name"].GetString());
            std::string typeName(jsonType["type"].GetString());
            assert(typeName == "donor" || typeName == "acceptor");
            auto probability = jsonType["probability"].GetDouble();
            HydrogenBondType bondType = typeName == "acceptor" ? HydrogenBondType::Acceptor : HydrogenBondType::Donor;
            HydrogenBondGeometry geometry = HydrogenBondGeometry::None;
            if (bondType == HydrogenBondType::Acceptor) {
                std::string geometryName(jsonType["geometry"].GetString());
                if (geometryName == "none") {
                    geometry = HydrogenBondGeometry::None;
                } else if (geometryName == "dir") {
                    geometry = HydrogenBondGeometry::Dir;
                } else if (geometryName == "cone") {
                    geometry = HydrogenBondGeometry::Cone;
                } else if (geometryName == "plane") {
                    geometry = HydrogenBondGeometry::Plane;
                } else {
                    REPORT(Reporter::FATAL) << "Unknown hydrogen bond geometry " << geometryName;
                    throw std::runtime_error("Unknown hydrogen bond geometry");
                }
            }
            std::string smarts(jsonType["smarts"].GetString());
            auto bondingType = std::make_shared<HydrogenBondingType>(bondType, name, probability, geometry, smarts);
            if (jsonType.HasMember("weight")) {
                auto weight = jsonType["weight"].GetDouble();
                bondingType->weight = weight;
            }
            hydrogenBondingTypes.push_back(bondingType);
        }

        auto const& jsonPartialCharges = d["partialChargeGroups"];
        partialCharges.clear();
        partialCharges.reserve(jsonPartialCharges.Size());
        for (auto& jsonRule: jsonPartialCharges.GetArray()) {
            std::string name(jsonRule["name"].GetString());
            std::string smarts(jsonRule["smarts"].GetString());
            auto formalCharge = jsonRule["formal"].GetInt();
            auto partialCharge = jsonRule["partial"].GetDouble();
            auto multiple = jsonRule["multiple"].GetBool();
            auto partialChargeGrp = std::make_shared<
                PartialCharge>(name, formalCharge, partialCharge, smarts, multiple);
            partialCharges.push_back(partialChargeGrp);
        }
    }
}
