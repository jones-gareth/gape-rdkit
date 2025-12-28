// Created by Gareth Jonews on 11/10/2024.
//


#pragma once

#include <set>
#include <Geometry/point.h>
#include <mol/Feature.h>

#include "SuperpositionChromosome.h"

namespace Gape {
    class FeaturePointSet;

    /**
     * A class to encode pharmacophore point scoring without reference to the base
     * molecule. Uses 3D relocation clustering of feature fitting points. When
     * scoring a cluster of fitting points each feature is tried in turn and the
     * best scoring feature is selected as a "base" feature.
     *
     * @author Gareth Jones
     *
     */
    class FeatureOverlay {
        const SuperpositionChromosome &superpositionChromosome;

        double groupFeatures(const FeatureType featureType);

        void setupFeaturePointSet(const FeatureType featureType);

        std::map<FeatureType, std::shared_ptr<FeaturePointSet> > featurePointSets;
        double donorHydrogenScore = 0.0;
        double acceptorAtomScore = 0.0;
        double aromaticRingScore = 0.0;
        double score;

    public:
        explicit FeatureOverlay(const SuperpositionChromosome &superpositionChromosome);

        int numberMolecules() const {
            return superpositionChromosome.superposition.numberMolecules();
        }

        const Superposition &getSuperposition() const { return superpositionChromosome.superposition; }

        double scoreOverlay();

        double getScore() const { return score; }

        double getDonorHydrogenScore() const { return donorHydrogenScore; }

        double getAcceptorAtomScore() const { return acceptorAtomScore; }

        double getAromaticRingScore() const { return aromaticRingScore; }

        const std::map<FeatureType, std::shared_ptr<FeaturePointSet> > &getFeaturePointSets() const {
            return featurePointSets;
        };

        int numberPharmacophorePoints() const;
    };

    /**
     * This class represents a cluster of pharmacophore features.
     *
     */
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

    /**
     * This class is used to store a set of feature points for a given feature
     * type
     *
     */
    class FeaturePointSet {
        const FeatureOverlay &featureOverlay;
        std::vector<std::shared_ptr<FeatureInformation> > features;
        std::set<std::shared_ptr<FeaturePoint> > featurePoints;
        std::set<std::shared_ptr<FeaturePoint> > freeFeaturePoints;
        int numberPharmacophorePoints = 0;
        double score = .0;
        const static int maxRelocations;


        int getPharmacophoreCount() const;

        void addFeature(FeatureInformation &feature, const double maxSqrDistance);

        bool relocate();

    public:
        const FeatureType featureType;

        FeaturePointSet(const FeatureOverlay &featureOverlay,
                        const std::vector<std::shared_ptr<FeatureInformation> > &features,
                        const FeatureType featureType): featureOverlay(featureOverlay),
                                                        features(features), featureType(featureType) {
            for ([[maybe_unused]] auto feature: features) {
                freeFeaturePoints.insert(std::make_shared<FeaturePoint>(featureOverlay));
            }
        }

        void groupPoints();

        double calculateScore();

        const std::vector<std::shared_ptr<FeatureInformation> > &getFeatures() const { return features; };

        friend class FeatureOverlay;
    };
}
