//
// Created by jones on 12/27/2022.
//

#ifndef GAPE_HYDROGENBONDINGTYPE_H
#define GAPE_HYDROGENBONDINGTYPE_H

#include <string>
#include <vector>

#include <GraphMol/GraphMol.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/SmilesParse/SmilesParse.h>

namespace Gape {

    using namespace RDKit;

    enum HydrogenBondType {
        Donor, Acceptor
    };

    enum HydrogenBondGeometry {
        None, Dir, Cone, Plane
    };

    struct HydrogenBondingType {
        const HydrogenBondType hydrogenBondType;
        const std::string name;
        const double probability;
        const HydrogenBondGeometry geometry;
        const std::string smarts;
        ROMol *query;
        double weight;

        HydrogenBondingType(const HydrogenBondType hydrogenBondType, const std::string &name, const double probability,
                            const HydrogenBondGeometry &geometry, const std::string &smarts);
        ~HydrogenBondingType();
    };

    using HydrogenBondingTypesList = std::vector<std::shared_ptr<const HydrogenBondingType>>;

    std::map<Atom *, std::shared_ptr<const HydrogenBondingType>>
    findHydrogenBondDonorsAndAcceptors(const HydrogenBondingTypesList &hydrogenBondingTypes, ROMol &mol);

} // Gape

#endif //GAPE_HYDROGENBONDINGTYPE_H
