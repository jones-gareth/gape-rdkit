//
// Created by jones on 11/10/2024.
//

#include "FeatureOverlay.h"

#include <util/Reporter.h>

#include "mol/MolUtil.h"

namespace Gape {
    void FeaturePoint::addFeaturePoint(const FeatureInformation& f, const bool updateCenter) {
        assert( features.find(&f) != features.end());
        features.insert(&f);
        const auto & point = f.point;
        sum.x += point.x;
        sum.y += point.y;
        sum.z += point.z;
        assert(checkSingleMolecules());
        if (updateCenter) {
            recalculateCenter();
        }
    }

    bool FeaturePoint::checkSingleMolecules() const {
        for (const auto & feature : features) {
            for (const auto &otherFeature: features) {
                if (feature == otherFeature) {
                    continue;
                }
                if (feature->feature->getMolecule() == otherFeature->feature->getMolecule()) {
                    return false;
                }
            }
        }
        return true;
    }

    void FeaturePoint::recalculateCenter() {
        const double numberPoints = features.size();
        center = sum/numberPoints;
    }

    void FeaturePoint::removeFeaturePoint(FeatureInformation& f, bool updateCenter) {
        assert( features.find(&f) != features.end());
        features.erase(&f);
        const auto & point = f.point;
        sum.x -= point.x;
        sum.y -= point.y;
        sum.z -= point.z;
        if (updateCenter) {
            recalculateCenter();
        }
    }

    double FeaturePoint::squareDistance(const FeatureInformation &otherFeatures) const {
        const auto & otherPoint = otherFeatures.point;
        double sqrDistance = squareDistance(point, otherPoint);
        return sqrDistance;
    }

    bool FeaturePoint::constainsMoleculeFeature(const SuperpositionMolecule &molecule) const {
        for (const auto & feature : features) {
            if (feature->feature->getMolecule() == &molecule) {
                return true;
            }
        }
        return false;
    }

    double FeaturePoint::scorePoint() {
        double maxScore = -std::numeric_limits<double>::max();
        baseFeature = nullptr;
        pharmacoporePoint = false;
        numberMatches = 0;
        score = .0;

        REPORT(Reporter::DEBUG) << "Entering FeaturePoint::score()";
        for (const auto & otherFeature : features) {
            double test = scoreFeature(*otherFeature);
            REPORT(Reporter::DEBUG) << "Test score: " << test;
            if (test > maxScore) {
                maxScore = test;
             REPORT(Reporter::DEBUG) << "Setting base feature to " << otherFeature->feature->;
            }
        }

    }





}
