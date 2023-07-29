#include "PartialCharge.h"
#include "util/Reporter.h"
#include <GraphMol/Substruct/SubstructMatch.h>

namespace Gape {
    const std::string PartialCharge::atomPartialChargeLabel = "__GapePartialCharge";
    const std::string PartialCharge::atomFormalChargeLabel = "__GapeFormalCharge";

    PartialCharge::PartialCharge(std::string n, const int fC, const double pC, std::string sma, const bool m) :
            name(std::move(n)), formalCharge(fC), partialCharge(pC), smarts(std::move(sma)), multiple(m) {
        try {
            query = SmartsToMol(smarts);
        }
        catch (...) {
            REPORT(Reporter::WARN) << "Failed to parse partial charge " << name << " pattern " << smarts;
        }
        if (query == nullptr) {
            REPORT(Reporter::WARN) << "Failed to parse partial charge " << name << " pattern " << smarts;
        }
    }

    void findPartialCharges(const PartialChargesList &partialChargesList, ROMol &mol) {
        SubstructMatchParameters params;
        params.uniquify = false;
        for (const auto &partialCharge: partialChargesList) {
            if (partialCharge->query != nullptr) {
                auto formalChargeSet = false;
                for (const auto &match: SubstructMatch(mol, *partialCharge->query, params)) {
                    const auto targetAtom = mol.getAtomWithIdx(match.front().second);
                    if (!targetAtom->hasProp(PartialCharge::atomPartialChargeLabel)) {
                        targetAtom->setProp(PartialCharge::atomPartialChargeLabel, partialCharge->partialCharge);
                        REPORT(Reporter::DETAIL) << "Setting partial charge on atom " << targetAtom->getIdx() << ":" <<
                                                 targetAtom->getSymbol() << " to " << partialCharge->partialCharge
                                                 << " using " << partialCharge->name << " definition";
                        if (partialCharge->multiple || !formalChargeSet) {
                            targetAtom->setProp(PartialCharge::atomFormalChargeLabel, partialCharge->formalCharge);
                            formalChargeSet = true;
                            REPORT(Reporter::DETAIL) << "Setting formal charge on atom " << targetAtom->getIdx() << ":"
                                                     <<
                                                     targetAtom->getSymbol() << " to " << partialCharge->formalCharge
                                                     << " using " << partialCharge->name << " definition";
                        }
                    }
                }
            }
        }
    }
} // namespace Gape
