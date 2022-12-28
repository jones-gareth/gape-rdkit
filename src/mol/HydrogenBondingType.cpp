//
// Created by jones on 12/27/2022.
//

#include "HydrogenBondingType.h"
#include "util/Reporter.h"
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/SmilesParse/SmilesParse.h>

namespace Gape {
    using namespace RDKit;

    HydrogenBondingType::HydrogenBondingType(const HydrogenBondType hydrogenBondType, const std::string &name,
                                             const double probability, const HydrogenBondGeometry &geometry,
                                             const std::string &smarts) : hydrogenBondType(hydrogenBondType),
                                                                          name(name), probability(probability),
                                                                          geometry(geometry), smarts(smarts) {
        try {
            query = SmartsToMol(this->smarts);
            weight = query->getNumAtoms();
        } catch (...) {
            REPORT(Reporter::WARN) << "Failed to parse hydrogen bond smarts pattern " << this->smarts;
        }
    }

    HydrogenBondingType::~HydrogenBondingType() {
        if (query != nullptr) {
            delete query;
        }
    }

    std::map<Atom *, std::shared_ptr<const HydrogenBondingType>>
    findHydrogenBondDonorsAndAcceptors(const HydrogenBondingTypesList &hydrogenBondingTypes, ROMol &mol) {
        std::map<Atom *, std::shared_ptr<const HydrogenBondingType>> features;
        SubstructMatchParameters params;
        params.uniquify = false;
        for (const auto &bondingType: hydrogenBondingTypes) {
            if (bondingType->query != nullptr) {
                for (const auto &match: SubstructMatch(mol, *bondingType->query, params)) {
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
} // Gape