// Created by Gareth Jonews on 11/10/2024.
//


#pragma once

#include <set>
#include <Geometry/point.h>
#include <mol/Feature.h>

#include "SuperpositionChromosome.h"

namespace Gape {
    class FeaturePoint {
        RDGeom::Point3D center;
        RDGeom::Point3D sum;
        double score = 0;
        std::set<FeatureInformation *> features;

        Feature *baseFeature;
        int numberMatches = 0;
        int testMatched = 0;
        bool pharmacoporePoint = false;

        bool checkSingleMolecules() const;
        double scoreFeature(FeatureInformation &otherFeature);
    public:
        void addFeaturePoint(const FeatureInformation& f, bool updateCenter);
        void recalculateCenter();
        void removeFeaturePoint(FeatureInformation& f, bool updateCenter);
        double squareDistance(const FeatureInformation& otherFeatures) const;
        bool constainsMoleculeFeature(const SuperpositionMolecule &molecule) const;
        double scorePoint();
    };

    class FeatureOverlay {
    };
}
