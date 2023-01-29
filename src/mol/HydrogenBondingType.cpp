//
// Created by jones on 12/27/2022.
//

#include "HydrogenBondingType.h"
#include "../util/Reporter.h"
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/SmilesParse/SmilesParse.h>

namespace Gape {
    using namespace RDKit;

    HydrogenBondingType::HydrogenBondingType(const HydrogenBondType hydrogenBondType, const std::string &name,
                                             const double probability, const HydrogenBondGeometry &geometry,
                                             const std::string &smartsIn) : hydrogenBondType(hydrogenBondType),
                                                                            name(name), probability(probability),
                                                                            geometry(geometry), smarts(smartsIn) {
        try {
            int wildcardNum = 0;
            for (int i = 0; i < smarts.length(); i++) {
                if (smarts[i] == '*') {
                    wildcardNum++;
                } else if (smarts[i] == 'Y') {
                    yWildcards.push_back(wildcardNum);
                    smarts[i] = '*';
                    wildcardNum++;
                } else if (smarts[i] == 'Z') {
                    zWildcards.push_back(wildcardNum);
                    smarts[i] = '*';
                    wildcardNum++;
                }
            }
            query = SmartsToMol(smarts);
            weight = query->getNumAtoms();
        } catch (...) {
            REPORT(Reporter::WARN) << "Failed to parse hydrogen bond smarts pattern " << this->smarts;
        }
    }

    bool HydrogenBondingType::matchZandY(const std::vector<std::pair<int, int>> &match, const ROMol &mol) const {
        if (yWildcards.empty() && zWildcards.empty()) {
            return true;
        }
        int wildCardNum = 0;
        int yAtomicWt = -1;
        int zAtomicWt = -1;
        for (const auto [queryIdx, molIdx]: match) {
            if (query->getAtomWithIdx(queryIdx)->getAtomicNum() == 0) {
                const auto molAtomicWt = mol.getAtomWithIdx(molIdx)->getAtomicNum();
                if (const auto yPresent = std::find(yWildcards.begin(), yWildcards.end(),
                                                    wildCardNum); yPresent != yWildcards.end()) {
                    if (yAtomicWt == -1) {
                        yAtomicWt = molAtomicWt;
                    } else if (yAtomicWt != molAtomicWt) {
                        return false;
                    }
                }
                if (const auto zPresent = std::find(zWildcards.begin(), zWildcards.end(),
                                                    wildCardNum); zPresent != zWildcards.end()) {
                    if (zAtomicWt == -1) {
                        zAtomicWt = molAtomicWt;
                    } else if (zAtomicWt != molAtomicWt) {
                        return false;
                    }
                }
                wildCardNum++;
            }
        }
        if (!yWildcards.empty()) {
            assert(yAtomicWt > 0);
        }
        if (!zWildcards.empty()) {
            assert(zAtomicWt > 0);
        }
        if (yAtomicWt > 0 && zAtomicWt > 0 && yAtomicWt == zAtomicWt) {
            return false;
        }
        return true;
    }

    HydrogenBondingType::~HydrogenBondingType() {
        if (query != nullptr) {
            delete query;
        }
    }

    std::map<Atom *, std::shared_ptr<const HydrogenBondingType>>
    findHydrogenBondDonorsOrAcceptors(const HydrogenBondType bondType, const HydrogenBondingTypesList &hydrogenBondingTypes, ROMol &mol) {
        std::map<Atom *, std::shared_ptr<const HydrogenBondingType>> features;
        SubstructMatchParameters params;
        params.uniquify = false;
        for (const auto &bondingType: hydrogenBondingTypes) {
            if (bondingType->hydrogenBondType != bondType) {
                continue;
            }
            if (bondingType->query != nullptr) {
                for (const auto &match: SubstructMatch(mol, *bondingType->query, params)) {
                    // wildcards can only match heavy atoms:
                    const auto hydrogenDummyMatch = std::find_if(match.begin(), match.end(),
                                                                 [&bondingType, &mol](const auto &pair) {
                                                                     auto queryAtom = bondingType->query->getAtomWithIdx(
                                                                             pair.first);
                                                                     auto molAtom = mol.getAtomWithIdx(pair.second);
                                                                     if ( molAtom->getAtomicNum() == 1) {
                                                                         assert(queryAtom->getAtomicNum() == 0);
                                                                         return true;
                                                                     }
                                                                     return false;
                                                                 });
                    if (hydrogenDummyMatch != match.end()) {
                        continue;
                    }
                    if (!bondingType->matchZandY(match, mol)) {
                        continue;
                    }
                    auto firstAtom = mol.getAtomWithIdx(match.front().second);
                    if (features.count(firstAtom)) {
                        if (features[firstAtom]->weight < bondingType->weight) {
                            features[firstAtom] = bondingType;
                        }
                    } else {
                        features[firstAtom] = bondingType;
                    }
                }
            }
        }

        for (const auto &[atom, bondingType]: features) {
            std::string typeStr(bondingType->hydrogenBondType == HydrogenBondType::Donor ? "Donor" : "Acceptor");
            REPORT(Reporter::DETAIL) << "Atom " << atom->getIdx() << " matches " << typeStr
                                     << " definition for " << bondingType->name;
        }

        return features;
    }


    std::map<Atom *, std::shared_ptr<const HydrogenBondingType>>
    findHydrogenBondDonors(const HydrogenBondingTypesList &hydrogenBondingTypes, ROMol &mol) {
        return findHydrogenBondDonorsOrAcceptors(HydrogenBondType::Donor, hydrogenBondingTypes, mol);
    }

    std::map<Atom *, std::shared_ptr<const HydrogenBondingType>>
    findHydrogenBondAcceptors(const HydrogenBondingTypesList &hydrogenBondingTypes, ROMol &mol) {
        return findHydrogenBondDonorsOrAcceptors(HydrogenBondType::Acceptor, hydrogenBondingTypes, mol);
    }
} // Gape