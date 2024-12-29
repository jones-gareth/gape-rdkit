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
        const FeatureOverlay &featureOverlay;

        /**
         * If this feature (in another molecule) is mappped to a feature in the base
         * molecule, we can store the base molecule feature here.
         */
        FeatureInformation *baseFeature = nullptr;
        int numberMatched = 0;

        /**
         * Set this true if we consider this feature to be part of a pharmacophore.
         */
        bool pharmacoporePoint = false;

        double scoreFeature(const FeatureInformation &testFeature, int &testMatched) const;

    public:
        explicit FeaturePoint(const FeatureOverlay &featureOverlay): featureOverlay(featureOverlay) {
        }

        bool checkSingleMolecules() const;

        void addFeaturePoint(FeatureInformation &f, bool updateCenter);

        void recalculateCenter();

        void removeFeaturePoint(FeatureInformation &f, bool updateCenter);

        double squareDistance(const FeatureInformation &otherFeature) const;

        bool containsMoleculeFeature(const SuperpositionMolecule &molecule) const;

        double scorePoint();

        void seed(FeatureInformation &feature);

        bool hasFeature(const FeatureInformation &feature) const;

        bool isPharmacophorePoint() const { return pharmacoporePoint; }

        /**
         * @return the number of features in this point.
         */
        int numberFeatures() const {
            return features.size();
        }
    };

    class FeaturePointSet {
        const FeatureOverlay &featureOverlay;
        const FeatureType featureType;
        std::set<FeatureInformation *> features;
        std::set<std::shared_ptr<FeaturePoint> > featurePoints;
        std::set<std::shared_ptr<FeaturePoint> > freeFeaturePoints;
        int numberPharmacophorePoints = 0;
        double score = .0;
        const static int maxRelocations;

        FeaturePointSet(const FeatureOverlay &featureOverlay, const FeatureType featureType,
                        const std::set<FeatureInformation *> &features): featureOverlay(featureOverlay),
                                                                         featureType(featureType), features(features) {
            for ([[maybe_unused]] auto feature: features) {
                freeFeaturePoints.insert(std::make_shared<FeaturePoint>(featureOverlay));
            }
        }

        double calculateScore();
        int getPharmacoporeCount() const;
        void groupPoints();
        void addFeature(FeatureInformation &feature, const double maxSqrDistance);
        bool relocate();
    };


    class FeatureOverlay {
        const SuperpositionGa &superpositionGa;

        

    public:
        explicit FeatureOverlay(const SuperpositionGa &superpositionGa): superpositionGa(superpositionGa) {
        }

        int numberMolecules() const {
            return superpositionGa.numberMolecules();
        }

        const SuperpositionGa &getSuperpositionGa() const { return superpositionGa; }
    };
}
