//
// Created by Gareth Jones on 12/27/2022.
//

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

    struct HydrogenBondingType;

    using HydrogenBondingTypesList = std::vector<std::shared_ptr<const HydrogenBondingType>>;

    struct HydrogenBondingType {
        const HydrogenBondType hydrogenBondType;
        const std::string name;
        const double probability;
        const HydrogenBondGeometry geometry;
        std::string smarts;
        ROMol *query;
        double weight;

        HydrogenBondingType(const HydrogenBondType hydrogenBondType, std::string name, const double probability,
                            const HydrogenBondGeometry &geometry, std::string smarts) ;

        ~HydrogenBondingType();

    private:
        std::vector<int> zWildcards, yWildcards;

        bool matchZandY(const std::vector<std::pair<int, int>> &match, const ROMol &mol) const;

        friend
        std::map<const Atom *, std::shared_ptr<const HydrogenBondingType>>
        findHydrogenBondDonorsOrAcceptors(const HydrogenBondType, const HydrogenBondingTypesList &hydrogenBondingTypes, ROMol &mol) ;
    };


    std::map<const Atom *, std::shared_ptr<const HydrogenBondingType>>
    findHydrogenBondDonors(const HydrogenBondingTypesList &hydrogenBondingTypes, ROMol &mol);

    std::map<const Atom *, std::shared_ptr<const HydrogenBondingType>>
    findHydrogenBondAcceptors(const HydrogenBondingTypesList &hydrogenBondingTypes, ROMol &mol);

} // Gape

