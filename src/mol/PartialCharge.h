#pragma once


#include <string>
#include <vector>

#include <GraphMol/GraphMol.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/SmilesParse/SmilesParse.h>

#include "../util/Reporter.h"

namespace Gape {

    using namespace RDKit;
    struct PartialCharge;

    using PartialChargesList = std::vector<std::shared_ptr<const PartialCharge>>;

    struct PartialCharge {
        const std::string name;
        const int formalCharge;
        const double partialCharge;
        std::string smarts;
        ROMol *query = nullptr;
        bool multiple;

        PartialCharge(std::string n, const int fC, const double pC, std::string sma, const bool m);

        friend void findPartialCharges(const PartialChargesList &partialChargesList, ROMol &mol);

        const static std::string atomPartialChargeLabel;
        const static std::string atomFormalChargeLabel;
    };

}
