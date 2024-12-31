//
// Created by Gareth Jones on 5/17/2024.
//

#pragma once

#include <vector>
#include <GraphMol/GraphMol.h>
#include "gape/SuperpositionMolecule.h"

using namespace RDKit;

namespace Gape {
    struct Gaussian {
        const double n;
        const double alpha;
        const RDGeom::Point3D point;

        Gaussian(const double n, const double alpha, const RDGeom::Point3D &point) : n(n), alpha(alpha), point(point) {}

        double volume() const;

        std::optional<Gaussian> intersection(const Gaussian &other) const;

        double overlapVolume(const Gaussian &other) const;
    };

    class GaussianList {
    private:
        std::vector<Gaussian> gaussians;

    public:
        GaussianList(const SuperpositionMolecule& molecule, const Conformer &conformer);

        GaussianList(const std::vector<Gaussian> &g) : gaussians(g) {}

        double overlapVolume(const Gaussian &g) const;

        double overlapVolume(const GaussianList &other) const;

        GaussianList intersection(const GaussianList &other) const;

        GaussianList intersection() const;

        double volume() const;

        // 2.7 is the the default atom volume
        static constexpr double atomicVolume = 2.7;

    };

}


