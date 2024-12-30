//
// Created by jones on 11/10/2024.
//

#include "FeatureOverlay.h"

#include <util/Reporter.h>

#include "mol/MolUtil.h"

namespace Gape {
    void FeaturePoint::addFeaturePoint(FeatureInformation &f, const bool updateCenter) {
        assert(features.find(&f) != features.end());
        features.insert(&f);
        const auto &point = f.point;
        sum.x += point.x;
        sum.y += point.y;
        sum.z += point.z;
        assert(checkSingleMolecules());
        if (updateCenter) {
            recalculateCenter();
        }
    }

    bool FeaturePoint::checkSingleMolecules() const {
        for (const auto &feature: features) {
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
        const double numberPoints = static_cast<int>(features.size());
        center = sum / numberPoints;
    }

    void FeaturePoint::removeFeaturePoint(FeatureInformation &f, bool updateCenter) {
        assert(features.find(&f) != features.end());
        features.erase(&f);
        const auto &point = f.point;
        sum.x -= point.x;
        sum.y -= point.y;
        sum.z -= point.z;
        if (updateCenter) {
            recalculateCenter();
        }
    }

    double FeaturePoint::squareDistance(const FeatureInformation &otherFeature) const {
        const auto &otherPoint = otherFeature.point;
        double sqrDistance = Gape::squareDistance(center, otherPoint);
        return sqrDistance;
    }

    bool FeaturePoint::containsMoleculeFeature(const SuperpositionMolecule &molecule) const {
        for (const auto &feature: features) {
            if (feature->feature->getMolecule() == &molecule) {
                return true;
            }
        }
        return false;
    }

    double FeaturePoint::scorePoint() {
        auto numMolecules = featureOverlay.numberMolecules();
        double maxScore = -std::numeric_limits<double>::max();
        baseFeature = nullptr;
        pharmacoporePoint = false;
        numberMatched = 0;
        score = .0;

        REPORT(Reporter::DEBUG) << "Entering FeaturePoint::score()";
        for (auto otherFeature: features) {
            int testMatched = 0;
            double test = scoreFeature(*otherFeature, testMatched);
            REPORT(Reporter::DEBUG) << "Test score: " << test;
            if (test > maxScore) {
                REPORT(Reporter::DEBUG) << "Setting base feature to " << otherFeature->feature->info();
                maxScore = test;
                baseFeature = otherFeature;
                numberMatched = testMatched;
            }

            // if we have two features the scores will be identical. Should
            // probably pick the feature with the highest group
            // probability.
            if (features.size() == 2)
                break;
        }

        assert(baseFeature != nullptr);

        // check for pharmacophore point
        if (numberMatched > 1 || (numMolecules == 2 && numberMatched == 1)) {
            assert(!baseFeature->isPharmacophoreFeature);
            baseFeature->isPharmacophoreFeature = true;
            baseFeature->numberMatches = numberMatched;
            pharmacoporePoint = true;
        }

        // Normalize to number of pairs
        maxScore = maxScore / (static_cast<double>(numMolecules) - 1.0);
        score = maxScore;
        return score;
    }

    /**
         * Test a feature as the base feature.
         *
         * @param testFeature
         * @return score that we get from using this feature as the base feature.
         */
    double FeaturePoint::scoreFeature(const FeatureInformation &testFeature, int &testMatched) const {
        double testScore = .0;
        testMatched = 0;

        auto testMolecule = testFeature.feature->getMolecule();
        const auto &settings = featureOverlay.getSuperpositionGa().getSuperposition().settings.getGapeParameters();

        // sum pair-wise score with all other features
        for (auto otherFeature: features) {
            if (otherFeature == &testFeature)
                continue;
            auto otherMolecule = otherFeature->feature->getMolecule();

            // only one feature per molecule can contribute to a
            // pharmacophore point.
            if (otherMolecule == testMolecule)
                continue;

            // get molecule weighting
            double weight = (testMolecule->getWeight() + otherMolecule->getWeight()) / 2.0;
            // determine feature pair score.
            const auto featureScore = testFeature.feature->score(*otherFeature->feature, testFeature.coordinates,
                                                                 otherFeature->coordinates);
            testScore += weight * featureScore.score;

            // count any "pharmacophore" matches
            if (settings.scalePharmacophore
                && featureScore.geometricScore > settings.geometricWeightCriterion)
                testMatched++;
        }

        // determine any pharmacophore scale up.
        if (testMatched > 2) {
            REPORT(Reporter::DEBUG) << "testMatched (2) " << testMatched;
            testScore *= pow(testMatched, settings.pharmacophoreFactor);
        }

        return testScore;
    }

    /**
     * Removes previous feature and scoring information from this point.
     * Seeds point with single feature. This should be the entry point for
     * creating new points.
     *
     * @param feature
     */
    void FeaturePoint::seed(FeatureInformation &feature) {
        features.clear();
        sum.x = sum.y = sum.z = 0;
        center.x = center.y = center.z = 0;
        numberMatched = 0;
        score = 0;
        pharmacoporePoint = false;
        addFeaturePoint(feature, true);
    }

    /**
  * @param feature
  * @return ture if this point contains the feature.
  */
    bool FeaturePoint::hasFeature(const FeatureInformation &feature) const {
        return features.find(const_cast<FeatureInformation *>(&feature)) != features.end();
    }

    double FeaturePointSet::calculateScore() {
        assert(featurePoints.size() > 0);
        assert(numberPharmacophorePoints == 0);


        // remove pharm point information here so that we always know
        // molecule pharm counts and feature counts are the same. This is
        // also done in groupPoints
        numberPharmacophorePoints = 0;
        for (const auto feature: features) {
            feature->isPharmacophoreFeature = false;
            feature->numberMatches = 0;
        }

        score = 0;
        for (auto featurePoint: featurePoints) {
            score += featurePoint->scorePoint();
            if (featurePoint->isPharmacophorePoint())
                numberPharmacophorePoints++;
        }
        return score;
    }

    int FeaturePointSet::getPharmacoporeCount() const {
        int count = 0;
        for (const auto feature: features) {
            if (feature->isPharmacophoreFeature) count++;
        }
        return count;
    }

    const int FeaturePointSet::maxRelocations = 100;

    /**
     * Groups all the points. First performs initial assignment of feature
     * points then does relocation until convergence.
     *
     * @see #relocate()
     * @see #addFeature(Feature)
     */
    void FeaturePointSet::groupPoints() {
        // clear out any previous grouping and scoring
        for (auto featurePoint: featurePoints) {
            freeFeaturePoints.insert(featurePoint);
        }
        featurePoints.clear();
        assert(featurePoints.size() + freeFeaturePoints.size() == features.size());

        // feature coordinates need to be this close (maxDistance) to form a
        // cluster
        auto maxDistance = Feature::getRadius() * 2;

        // features from the same
        // molecule will always occupy their own cluster
        double maxSqrDistance = maxDistance * maxDistance;

        // / initial assignment
        for (auto feature: features)
            addFeature(*feature, maxSqrDistance);

        // relocations
        int nRelocations = 0;
        while (relocate()) {
            nRelocations++;
            if (nRelocations == maxRelocations) {
                REPORT(Reporter::INFO) << "groupPoints no convergence after " << maxRelocations << " relocations";
                break;
            }
        }

        for (auto featurePoint: featurePoints)
            assert(featurePoint->checkSingleMolecules());
    }

    bool FeaturePointSet::relocate() {
        bool relocate = false;

        assert(featurePoints.size() + freeFeaturePoints.size() == features.size());

        for (auto feature: features) {
            // foreach feature find the closest featurePoint and the point
            // containing that feature.
            std::shared_ptr<FeaturePoint> closestPoint = nullptr;
            std::shared_ptr<FeaturePoint> currentPoint = nullptr;
            double minDistance = std::numeric_limits<double>::max();
            for (auto featurePoint: featurePoints) {
                double distance = featurePoint->squareDistance(*feature);
                if (distance < minDistance) {
                    minDistance = distance;
                    closestPoint = featurePoint;
                }
                if (featurePoint->hasFeature(*feature))
                    currentPoint = featurePoint;
            }

            assert(closestPoint != nullptr);
            assert(currentPoint != nullptr);

            // check if feature is already the closest point.
            if (closestPoint == currentPoint)
                continue;

            // check constraint that each point only contains one feature
            // from a molecule
            if (closestPoint->containsMoleculeFeature(*feature->feature->getMolecule()))
                continue;

            // otherwise move it to the closest point.
            relocate = true;
            currentPoint->removeFeaturePoint(*feature, false);
            closestPoint->addFeaturePoint(*feature, false);

            // if the relocation leaves an empty point remove it.
            if (currentPoint->numberFeatures() == 0) {
                featurePoints.erase(currentPoint);
                freeFeaturePoints.insert(currentPoint);
                assert(featurePoints.size() + freeFeaturePoints.size() == features .size());
            }
        }

        // update centers
        for (auto featurePoint: featurePoints)
            featurePoint->recalculateCenter();

        return relocate;
    }

    /**
     * Adds a feature to the initial grouping.
     *
     * @param feature
     */
    void FeaturePointSet::addFeature(FeatureInformation &feature, const double maxSqrDistance) {
        const auto molecule = feature.feature->getMolecule();
        bool added = false;

        // loop though any existing points
        for (auto featurePoint: featurePoints) {
            double sqrDistance = featurePoint->squareDistance(feature);
            if (featurePoint->containsMoleculeFeature(*molecule))
                continue;
            // criteria if the current point does not contain any
            // feature from this molecule.
            if (sqrDistance < maxSqrDistance) {
                // add updating center
                featurePoint->addFeaturePoint(feature, true);
                added = true;
                break;
            }
        }

        if (!added) {
            // no close points so create a new feature point with this
            // feature.
            auto newFeaturePoint = *std::next(freeFeaturePoints.begin());
            freeFeaturePoints.erase(newFeaturePoint);
            newFeaturePoint->seed(feature);
            featurePoints.insert(newFeaturePoint);
        }
    }

    /**
     * Creates an object for doing the feature overlay. Initializes
     * FeaturePointSets for built-in and user features. Sets the pharmacophore
     * counts for molecules to zero. If you're re-scoring (such as in a GA run)
     * you need to set the molecule pharmacophore counts to zero before calling
     * scoring functions.
     *
     * @param problem
     */
    FeatureOverlay::FeatureOverlay(const SuperpositionChromosome &superpositionChromosome): superpositionChromosome(superpositionChromosome) {
        assert(numberMolecules() > 2);
        const auto& settings = superpositionChromosome.superpositionGa.getSuperposition().settings.getGapeParameters();

        if (settings.donorHydrogenWeight > 0) {
            setup(FeatureType::DonorInteractionPoint);
        }
        if (settings.acceptorAtomWeight > 0) {
            setup(FeatureType::AcceptorAtomFeature);
        }
        if (settings.aromaticRingWeight > 0) {
            setup(FeatureType::AromaticRing);
        }
    }

	/**
	 * Initializes a feature set. Creates a FeaturePointSet.
	 *
	 * @param featureType
	 */
    void FeatureOverlay::setup(const FeatureType featureType) {
        std::set<std::shared_ptr<FeatureInformation>> setFeatures;
        auto moleculeNumber = 0;
        const auto& superpositionCoordinates = superpositionChromosome.getFittedCoordinates();
        for (const auto& molecule: superpositionChromosome.superpositionGa.getSuperposition().getMolecules()) {
            const auto& coordinates = superpositionCoordinates[moleculeNumber];
            const auto& features = molecule->getFeatures().at(featureType);
            for (const auto& feature: features) {
                auto featureInformation= std::make_shared<FeatureInformation>(feature.get(), *coordinates);
                setFeatures.insert(featureInformation);
            }
            moleculeNumber++;
        }
        auto featurePointSet = std::make_shared<FeaturePointSet>(*this, featureType, setFeatures);
        featurePointSets.emplace(featureType, featurePointSet);
    }

	/**
	 * Takes the current saved feature coordinates and determines GAPE score.
	 * More typically we'll call scoreFeature for each feature.
	 *
	 * @return feature score.
	 */
    double FeatureOverlay::scoreOverlay() {
        score = .0;
        const auto& settings = superpositionChromosome.superpositionGa.getSuperposition().settings.getGapeParameters();

        if (settings.donorHydrogenWeight > 0) {
            score += settings.donorHydrogenWeight * groupFeatures(FeatureType::DonorInteractionPoint);
        }
        if (settings.acceptorAtomWeight > 0) {
            score += settings.acceptorAtomWeight * groupFeatures(FeatureType::AcceptorAtomFeature);
        }
        if (settings.aromaticRingWeight > 0) {
            score += settings.aromaticRingWeight * groupFeatures(FeatureType::AromaticRing);
        }

        return score;
    }

	/**
	 * Determines GAPE score for a feature. First clusters in 3D space then
	 * scores clusters.
	 *
	 * @param featureSetNo
	 * @return
	 */
    double FeatureOverlay::groupFeatures(const FeatureType featureType) {
        auto featurePointSet = featurePointSets[featureType];
        assert(featurePointSet->featureType == featureType);
        featurePointSet->groupPoints();
        const auto setScore = featurePointSet->calculateScore();
        return setScore;
    }

}
