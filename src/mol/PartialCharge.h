#pragma once


#include <string>
#include <vector>

#include <GraphMol/GraphMol.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/SmilesParse/SmilesParse.h>

namespace Gape {

    using namespace RDKit;

    enum class HydrogenBondType {
        Donor, Acceptor
    };

    enum class HydrogenBondGeometry {
        None, Dir, Cone, Plane
    };

    struct PartialCharge;

    using PartialChargesList = std::vector<std::shared_ptr<const PartialCharge>>;

    struct PartialCharge {
        const std::string name;
        const double formalCharge;
        const double partialCharge;
        const HydrogenBondGeometry geometry;
        std::string smarts;
        ROMol* query;
        bool multiple;

        


    }

}
