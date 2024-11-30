//
// Created by Gareth Jones on 5/17/2024.
//

#include "GaussianList.h"
#include "Reporter.h"
#include <GraphMol/PeriodicTable.h>

namespace Gape {

    const double VDW_CUTOFF = 1.0e-2;

    GaussianList::GaussianList(const Gape::SuperpositionMolecule molecule, const RDKit::Conformer &conformer) {
        gaussians.clear();
        const auto &hydrophobicFeatures = molecule.getFeatures().at(FeatureType::HydrophobicAtom);
        gaussians.reserve(hydrophobicFeatures.size());
        for (const auto &feature: hydrophobicFeatures) {
            const auto atom = feature->getAtom();
            const auto atomIdx = atom->getIdx();
            const auto &point = conformer.getAtomPos(atomIdx);


        }
    }

    /**
     * Determines the alpha parameter (or Gaussian width paramterer) for an
     * atomic Gaussian.
     *
     * @return
     */
    double gaussianWidth(int atomicNumber) {
        const auto *ptable = RDKit::PeriodicTable::getTable();
        auto r = ptable->getRvdw(atomicNumber);
        // 2.7 is the the default atom volume
        auto alpha = pow((3.0 * 2.7) / (4 * M_PI * r * r * r),
                         2.0 / 3.0) * M_PI;
        return alpha;
    }

    double gaussianVolume(const double n, const double alpha) {
        auto f = M_PI / alpha;
        auto vol = n * f * sqrt(f);
        return vol;
    }

    double Gaussian::volume() const {
        return gaussianVolume(n, alpha);
    }

    double Gaussian::overlapVolume(const Gaussian &other) const {
        auto rSqr = (point - other.point).lengthSq();
        auto alpha2 = alpha + other.alpha;
        auto val = -((alpha * other.alpha) / alpha2) * rSqr;
        auto n2 = n * other.n * exp(val);
        return gaussianVolume(n2, alpha2);
    }

    /**
     * Intersects two Gausians to create ane Gaussian in this list. The Gausian
     * is only created if there is a minimum overlap.
     */
    std::optional<Gaussian> Gaussian::intersection(const Gape::Gaussian &other) const {
        const auto rSqr = (point - other.point).lengthSq();
        double alpha3 = alpha + other.alpha;
        const auto val = -((alpha * other.alpha) / alpha3) * rSqr;
        double n3 = n * other.n * exp(val);
        double vol = gaussianVolume(n3, alpha3);
        // intersectVolume += vol;

        if (vol < VDW_CUTOFF)
            return std::nullopt;

        // new gaussian center
        const auto point1 = point * alpha;
        const auto point2 = other.point * other.alpha;
        const auto point3 = (point1 + point2) / (alpha * other.alpha);

        Gaussian rtn(n3, alpha3, point3);
        return rtn;
    }

    double GaussianList::overlapVolume(const Gape::Gaussian &g) const {
        double vol = 0;
        for (const auto &g2: gaussians) {
            vol += g2.overlapVolume(g);
        }
        REPORT(Reporter::DEBUG) << "volume " << vol;
        return vol;
    }

    /**
     * @param b
     * @return a new Gaussian list by intersecting this list with another
     *         Gaussian list. Each Gaussian in this list is intersected with
     *         every Gaussian in the other list. Used a a first step in
     *         determining intersection volumes.
     */
    GaussianList GaussianList::intersection(const Gape::GaussianList &other) const {
        std::vector<Gaussian> intersection;
        for (const auto &g1: gaussians) {
            for (const auto &g2: other.gaussians) {
                const auto g3 = g1.intersection(g2);
                if (g3.has_value()) {
                    intersection.push_back(g3.value());
                }
            }
        }

        REPORT(Reporter::DEBUG) << "Overlay intersection size " << intersection.size();
        GaussianList rtn(intersection);
        return rtn;
    }

    /**
     * @return a new Gaussian list by intersecting this list with itself. Each
     *         Gaussian is intersected with every other Gaussian, but not with
     *         itself. Use to determine corrections in volume calculations.
     */
    GaussianList GaussianList::intersection() const {
        const auto size = gaussians.size();
        std::vector<Gaussian> intersection;
        REPORT(Reporter::DEBUG) << "Max number of self intersections gaussians " << size * (size - 1) / 2;
        for (int i = 0; i < size; i++) {
            const auto g1 = gaussians[i];
            for (int j = i + 1; j < size; j++) {
                const auto g2 = gaussians[j];
                const auto g3 = g1.intersection(g2);
                if (g3.has_value()) {
                    intersection.push_back(g3.value());
                }
            }
        }
        REPORT(Reporter::DEBUG) << "Actual overlay intersection size " << intersection.size();
        GaussianList rtn(intersection);
        return rtn;
    }

    /**
     * @param b
     * @return The intersection volume (simple sum of volume of intersection
     *         Gaussians) determined by intersecting each Gaussian in this list
     *         is intersected with every Gaussian in the other list.
     */
    double GaussianList::overlapVolume(const Gape::GaussianList &other) const {
        auto vol = .0;
        for (const auto &g1: gaussians) {
            for (const auto &g2: other.gaussians) {
                vol += g1.overlapVolume(g2);
            }
        }
        REPORT(Reporter::DEBUG) << "Overlap of " << gaussians.size() << " with "
        << other.gaussians.size() << " volume " << vol;
        return vol;
    }

}