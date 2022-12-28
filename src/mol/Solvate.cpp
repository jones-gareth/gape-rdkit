//
// Created by Gareth Jones on 11/27/2022.
//

#include "Solvate.h"
#include "../util/Reporter.h"
#include <GraphMol/Substruct/SubstructMatch.h>

namespace Gape {
    using namespace RDKit;

    void solvateMolecule(const SolvationRuleList &rules, RWMol &mol) {
        // the molecule should have had hydrogens added prior to solvation

        // the original documentation for this process requires that
        // we do solvation twice- the second pass can apply any rule
        // corrections.  Not sure what that means, but have left the
        // two passes in;

        for (int pass = 0; pass < 2; ++pass) {
            REPORT(Reporter::DETAIL) << "Solvation pass " << pass + 1;
            bool matched = false;

            for (const auto &solvationRule: rules) {
                std::string typeStr(solvationRule->solvationType == SolvationType::Acid ? "acid" : "base");
                for (const auto &match: SubstructMatch(mol, *solvationRule->query)) {
                    REPORT(Reporter::DETAIL) << "Solvation rule " << solvationRule->name << " for " << typeStr
                                             << " matches molecule " << mol.getProp<std::string>("_Name");
                    matched = true;
                    auto firstAtom = mol.getAtomWithIdx(match.front().second);
                    switch (solvationRule->solvationType) {
                        case SolvationType::Acid: {
                            // remove neighbor hydrogen and set first matching atom -1 charge
                            Atom *hydrogen = nullptr;
                            for (auto &neighbor: mol.atomNeighbors(firstAtom)) {
                                if (neighbor->getAtomicNum() == 1) {
                                    hydrogen = neighbor;
                                    break;
                                }
                            }
                            if (hydrogen != nullptr) {
                                mol.removeAtom(hydrogen);
                                firstAtom->setFormalCharge(-1);
                                firstAtom->updatePropertyCache(false);
                            } else {
                                REPORT(Reporter::WARN) << "Failed to find hydrogen for matching acid group "
                                                       << solvationRule->name;
                            }
                        }
                            break;
                        case SolvationType::Base: {
                            // Add hydrogen and set first matching atom +1 charge
                            const auto newHydrogen = new Atom(1);
                            const auto hydrogenIndex = mol.addAtom(newHydrogen);
                            const auto bond = new Bond(Bond::BondType::SINGLE);
                            bond->setBeginAtomIdx(firstAtom->getIdx());
                            bond->setEndAtomIdx(hydrogenIndex);
                            mol.addBond(bond);
                            firstAtom->setFormalCharge(1);
                            firstAtom->updatePropertyCache(false);
                        }
                            break;
                    }
                }
            }

            if (!matched) {
                break;
            }
        }
    }
}