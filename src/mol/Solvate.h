//
// Created by Gareth Jones on 11/27/2022.
//

#ifndef GAPE_SOLVATE_H
#define GAPE_SOLVATE_H

#include <string>
#include <GraphMol/GraphMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>


namespace Gape {
    using namespace RDKit;

    enum SolvationType {
        Acid, Base
    };

    struct SolvationRule {
        const SolvationType solvationType;
        const std::string name;
        const std::string smarts;
        const double pKa;
        ROMol *query;

        SolvationRule(const SolvationType &solvationType, std::string name, std::string smarts, double pKa)
                :
                solvationType(solvationType), name(std::move(name)), smarts(std::move(smarts)), pKa(pKa) {
            query = SmartsToMol(this->smarts);
        }

        ~SolvationRule() {
            delete query;
        }
    };

    using SolvationRuleList = std::vector<std::shared_ptr<const SolvationRule>>;

    void solvateMolecule(const SolvationRuleList &rules, RWMol &mol);

}
#endif //GAPE_SOLVATE_H
