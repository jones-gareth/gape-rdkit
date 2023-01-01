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

    class HydrogenBondingType;

    using HydrogenBondingTypesList = std::vector<std::shared_ptr<const HydrogenBondingType>>;

    struct HydrogenBondingType {
        const HydrogenBondType hydrogenBondType;
        const std::string name;
        const double probability;
        const HydrogenBondGeometry geometry;
        std::string smarts;
        ROMol *query;
        double weight;

        HydrogenBondingType(const HydrogenBondType hydrogenBondType, const std::string &name, const double probability,
                            const HydrogenBondGeometry &geometry, const std::string &smarts) ;

        ~HydrogenBondingType();

    private:
        std::vector<int> zWildcards, yWildcards;

        bool matchZandY(const std::vector<std::pair<int, int>> &match, const ROMol &mol) const;

        friend
        std::map<Atom *, std::shared_ptr<const HydrogenBondingType>>
        findHydrogenBondDonorsOrAcceptors(const HydrogenBondType, const HydrogenBondingTypesList &hydrogenBondingTypes, ROMol &mol) ;
    };


    std::map<Atom *, std::shared_ptr<const HydrogenBondingType>>
    findHydrogenBondDonors(const HydrogenBondingTypesList &hydrogenBondingTypes, ROMol &mol);

    std::map<Atom *, std::shared_ptr<const HydrogenBondingType>>
    findHydrogenBondAcceptors(const HydrogenBondingTypesList &hydrogenBondingTypes, ROMol &mol);

} // Gape

#endif //GAPE_HYDROGENBONDINGTYPE_H
