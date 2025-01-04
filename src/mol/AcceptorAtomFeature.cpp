//
// Created by Gareth Jones on 3/23/2023.
//

#include <boost/format.hpp>
#include <boost/algorithm/string.hpp>
#include "AcceptorAtomFeature.h"

#include "HydrogenBondingType.h"
#include "util/Reporter.h"
#include "util/TransformOps.h"
#include "mol/PartialCharge.h"

namespace Gape {
    thread_local double AcceptorAtomFeature::maxLonePairLonePairAngle = 50.0 * M_PI / 180.0;
    thread_local double AcceptorAtomFeature::minLonePairLonePairAngle = 20.0 * M_PI / 180.0;
    thread_local double AcceptorAtomFeature::maxForwardAcceptorAngle = 80.0 * M_PI / 180.0;
    thread_local double AcceptorAtomFeature::minForwardAcceptorAngle = 50.0 * M_PI / 180.0;
    thread_local double AcceptorAtomFeature::maxPlanePlaneAngle = 50.0 * M_PI / 180.0;
    thread_local double AcceptorAtomFeature::minPlanePlaneAngle = 20.0 * M_PI / 180.0;
    thread_local double AcceptorAtomFeature::maxPlaneLonePairAngle = 50.0 * M_PI / 180.0;
    thread_local double AcceptorAtomFeature::minPlaneLonePairAngle = 20.0 * M_PI / 180.0;
    thread_local double AcceptorAtomFeature::hBondLen = 2.9;
    thread_local double AcceptorAtomFeature::chargeFactor = 2.0;
    thread_local double AcceptorAtomFeature::matchFactor = 1.0;
    double AcceptorAtomFeature::maxAcceptorPartialCharge = -0.2;

    namespace Detail {
        /**
         * Scoring function for an angle constraint. Returns 1 if angle < min, 0 if
         * angle > max and linear interpolates otherwise.
         *
         * @param angle
         * @param max
         * @param min
         * @return
         */
        double scoreAngle(double angle, double max, double min) {
            if (angle > max)
                return .0;
            if (angle < min)
                return 1.0;
            return 1 - (angle - min) / (max - min);
        }

        void getPlaneNormal(const RDGeom::Point3D &point, const std::vector<RDGeom::Point3D> &lonePairs,
                            RDGeom::Point3D &normal) {
            assert(lonePairs.size() == 2);
            const auto v1 = lonePairs[0] - point;
            const auto v2 = lonePairs[1] - point;
            normal = v1.crossProduct(v2);
        }

        /**
         * For two acceptors that each accept along a single lone pair, use angle
         * constraints to check that they have compatible orientation.
         */
        double lonePairLonePairScore(const RDGeom::Point3D &coordinate, const RDGeom::Point3D &otherCoordinate,
                                     const std::vector<RDGeom::Point3D> &lonePairs,
                                     const std::vector<RDGeom::Point3D> &otherLonePairs) {
            assert(lonePairs.size() == 1);
            assert(otherLonePairs.size() == 1);

            const auto &lonePair = lonePairs[0];
            const auto &otherLonePair = otherLonePairs[0];

            auto angle = angleBetween(lonePair, coordinate, otherLonePair, otherCoordinate);
            if (angle > M_PI)
                angle = M_PI;
            REPORT(Reporter::DEBUG) << "LP LP angle " << (angle * 180.0 / M_PI);
            const auto score = scoreAngle(angle, AcceptorAtomFeature::maxLonePairLonePairAngle,
                                          AcceptorAtomFeature::minLonePairLonePairAngle);
            return score;
        }


        /**
         * For two acceptors that each accept along a plane of two lone pairs, use
         * angle constraints to check that they have compatible orientation.
         */
        double planePlaneScore(const RDGeom::Point3D &coordinate, const RDGeom::Point3D &otherCoordinate,
                               const std::vector<RDGeom::Point3D> &lonePairs,
                               const std::vector<RDGeom::Point3D> &otherLonePairs) {
            assert(lonePairs.size() == 2);
            assert(otherLonePairs.size() == 2);

            RDGeom::Point3D normal, otherNormal;
            getPlaneNormal(coordinate, lonePairs, normal);
            getPlaneNormal(otherCoordinate, otherLonePairs, otherNormal);
            auto angle = normal.angleTo(otherNormal);
            if (angle > M_PI / 2.0) angle = M_PI - angle;
            REPORT(Reporter::DEBUG) << "Plane plane angle " << (angle * 180.0 / M_PI);
            const auto score = scoreAngle(angle, AcceptorAtomFeature::maxPlanePlaneAngle,
                                          AcceptorAtomFeature::minPlanePlaneAngle);
            return score;
        }

        double planeLonePairScore(const RDGeom::Point3D &coordinate, const RDGeom::Point3D &otherCoordinate,
                                  const std::vector<RDGeom::Point3D> &lonePairs,
                                  const std::vector<RDGeom::Point3D> &otherLonePairs) {
            assert(lonePairs.size() == 2);
            assert(otherLonePairs.size() == 1);

            RDGeom::Point3D normal;
            getPlaneNormal(coordinate, lonePairs, normal);
            const auto v = otherLonePairs[0] - otherCoordinate;
            auto angle = normal.angleTo(v);
            angle = angle - M_PI / 2.0;
            if (angle < 0)
                angle = -angle;
            const auto score = scoreAngle(angle, AcceptorAtomFeature::maxPlaneLonePairAngle,
                                          AcceptorAtomFeature::minPlaneLonePairAngle);
            return score;
        }
    }

    void AcceptorAtomFeature::getSolvationPoint(const SuperpositionCoordinates &superpositionCoordinates,
                                                RDGeom::Point3D &solvationPoint) const {
        const auto &lonePairs = superpositionCoordinates.getFeatureCoordinates(FeatureType::AcceptorAtomFeature, atom);
        if (lonePairs.size() == 1) {
            solvationPoint = lonePairs[0];
        } else if (lonePairs.size() == 2) {
            const auto midPoint = (lonePairs[0] + lonePairs[1]) / 2.0;
            const auto &coordinate = superpositionCoordinates.getConformer().getAtomPos(atom->getIdx());
            auto diff = midPoint - coordinate;
            diff.normalize();
            diff *= hBondLen;
            solvationPoint = coordinate + diff;
        } else {
            assert(false);
        }
    }


    std::vector<std::shared_ptr<Feature> >
    AcceptorAtomFeature::findAcceptorAtoms(const SuperpositionMolecule *superpositionMolecule) {
        auto &mol = superpositionMolecule->getMol();
        auto &acceptors = superpositionMolecule->getAcceptors();
        std::vector<std::shared_ptr<Feature> > features;

        int featureSetNumber = 0;
        for (const auto& [atom, t]: acceptors) {
            auto acceptorAtom = mol.getAtomWithIdx(atom->getIdx());
            auto feature = std::make_shared<AcceptorAtomFeature>(featureSetNumber++, superpositionMolecule,
                                                                 acceptorAtom);
            features.push_back(std::static_pointer_cast<Feature>(feature));
        }
        return features;
    }

    AcceptorAtomFeature::AcceptorAtomFeature(const int featureSetNum, const SuperpositionMolecule *spMol,
                                             const Atom *featureAtom): AcceptorAtomFeature(featureSetNum) {
        molecule = spMol;
        atom = featureAtom;
        const auto &acceptors = molecule->getAcceptors();
        const auto acceptorIt = acceptors.find(atom);
        assert(acceptorIt != acceptors.end());
        hydrogenBondingType = acceptorIt->second.get();
        charged = PartialCharge::getPartialCharge(featureAtom) < maxAcceptorPartialCharge;
        REPORT(Reporter::DEBUG) << "N Lone Pairs " << nLonePairs;
        acceptorAtom = std::make_unique<AcceptorAtom>(molecule, atom);
    }

    std::string AcceptorAtomFeature::pharmLabel() const {
        auto name = boost::to_upper_copy<std::string>(hydrogenBondingType->name);
        std::replace(name.begin(), name.end(), ' ', '_');
        std::string options;
        if (charged) {
            options += " [charged = yes]";
        }
        const auto label = (boost::format("ATOM_%d") % (atom->getIdx() + 1)).str();
        return featureSetName + " " + label + options;
    }

    std::string AcceptorAtomFeature::info() const {
        const auto format = boost::format("Acceptor [%d]") % (atom->getIdx() + 1);
        return format.str();
    }

    const RDGeom::Point3D &AcceptorAtomFeature::getFittingPoint(
        const SuperpositionCoordinates &superpositionCoordinates) const {
        return superpositionCoordinates.getConformer().getAtomPos(atom->getIdx());
    }

    void AcceptorAtomFeature::calculateCoordinates(SuperpositionCoordinates &superpositionCoordinates) const {
        std::vector<RDGeom::Point3D> lonePairs;
        const auto &conformer = superpositionCoordinates.getConformer();
        acceptorAtom->addLonePairs(conformer, lonePairs);
        superpositionCoordinates.addFeatureCoordinates(FeatureType::AcceptorAtomFeature, atom, lonePairs);
    }

    FeatureScore AcceptorAtomFeature::score(const Feature &otherFeature, const SuperpositionCoordinates &coordinates,
                                            const SuperpositionCoordinates &otherCoordinates) const {
        const auto &other = dynamic_cast<const AcceptorAtomFeature &>(otherFeature);

        // Score for acceptor overlay
        const auto &coordinate = coordinates.getConformer().getAtomPos(atom->getIdx());
        const auto &otherCoordinate = otherCoordinates.getConformer().getAtomPos(other.atom->getIdx());
        const auto sqrDistance = (coordinate - otherCoordinate).lengthSq();
        auto accScore = Feature::score(sqrDistance);
        if (accScore > maximumGaussianScore) {
            accScore = maximumGaussianScore;
        }

        // Corrections to make sure fitting point is solvent accessible
        double correct = .0;
        RDGeom::Point3D solvationPoint, otherSolvationPoint;
        getSolvationPoint(coordinates, solvationPoint);
        other.getSolvationPoint(otherCoordinates, otherSolvationPoint);
        auto midPoint = (solvationPoint + otherSolvationPoint) / 2.0;
        double molVol = solvationPenalty(midPoint, coordinates);
        double otherMolVol = other.solvationPenalty(midPoint, otherCoordinates);

        correct = molVol > otherMolVol ? molVol : otherMolVol;
        if (correct < .0)
            correct = .0;
        accScore -= correct;
        if (accScore <= 0) {
            FeatureScore result(.0, .0);
            return result;
        }

        // Both acceptors pointing in the same direction?
        const auto angle = angleBetween(solvationPoint, coordinate, otherSolvationPoint, otherCoordinate);
        const auto forwardScore = Detail::scoreAngle(angle, maxForwardAcceptorAngle, minForwardAcceptorAngle);
        if (forwardScore <= .0) {
            FeatureScore result(.0, .0);
            return result;
        }

        // Compatible geometries?
        double geometryScore = .0;
        const auto geometry = hydrogenBondingType->geometry;
        const auto otherGeometry = other.hydrogenBondingType->geometry;
        const auto &lonePairs = coordinates.getFeatureCoordinates(FeatureType::AcceptorAtomFeature, atom);
        const auto &otherLonePairs = otherCoordinates.
                getFeatureCoordinates(FeatureType::AcceptorAtomFeature, other.atom);

        if (geometry == HydrogenBondGeometry::None || otherGeometry == HydrogenBondGeometry::None || geometry ==
            HydrogenBondGeometry::Cone || otherGeometry == HydrogenBondGeometry::Cone)
            // For acceptors that accept in a cone or have no directionality
            // assume that the general forward constraint is sufficient.
            geometryScore = 1.0;

        else if (geometry == HydrogenBondGeometry::Dir && otherGeometry == HydrogenBondGeometry::Dir)
            geometryScore = Detail::lonePairLonePairScore(coordinate, otherCoordinate, lonePairs, otherLonePairs);

        else if (geometry == HydrogenBondGeometry::Plane && otherGeometry == HydrogenBondGeometry::Plane)
            geometryScore = Detail::planePlaneScore(coordinate, otherCoordinate, lonePairs, otherLonePairs);

        else if (geometry == HydrogenBondGeometry::Plane && otherGeometry == HydrogenBondGeometry::Dir)
            geometryScore = Detail::planeLonePairScore(coordinate, otherCoordinate, lonePairs, otherLonePairs);

        else if (geometry == HydrogenBondGeometry::Dir && otherGeometry == HydrogenBondGeometry::Plane)
            // ignore swapped argument warning
            geometryScore = Detail::planeLonePairScore(otherCoordinate, coordinate, otherLonePairs, lonePairs);

        assert(!isnan(geometryScore));

        double score = accScore * forwardScore * geometryScore;
        const auto geometricScore = score;
        // Type matching scale up
        double prob = hydrogenBondingType->probability + other.hydrogenBondingType->probability;
        score *= prob;
        if (hydrogenBondingType->name == other.hydrogenBondingType->name)
            score *= matchFactor;

        // charge scale up
        if (charged && other.charged)
            score *= chargeFactor;

        REPORT(Reporter::DEBUG) << info() << " " << other.info() << " score " << score << " geometric score " <<
			 geometricScore << " geometry score " << geometryScore << " forward score " << forwardScore <<
             " solvation correction " << correct << " accScore " << accScore;

        const FeatureScore result(score, geometricScore);
        return result;
    }

    std::unique_ptr<PharmFeatureGeometry> AcceptorAtomFeature::getPharmFeatureGeometry(
        const SuperpositionCoordinates &superpositionCoordinates) const {
        const auto &coord = superpositionCoordinates.getConformer().getAtomPos(atom->getIdx());
        const auto &lonePairs = superpositionCoordinates.getFeatureCoordinates(FeatureType::AcceptorAtomFeature, atom);

        switch (hydrogenBondingType->geometry) {
            case HydrogenBondGeometry::Dir:
                assert(lonePairs.size() == 1);
                return std::make_unique<VectorPharmFeatureGeometry>(coord, lonePairs[0]);
            case HydrogenBondGeometry::Plane:
                assert(lonePairs.size() == 2);
                return std::make_unique<ArcFeatureGeometry>(coord, lonePairs[0], lonePairs[1]);
            case HydrogenBondGeometry::Cone: {
                // Cone acceptor is defined by a single lone pair , but for the cone
                // feature we want two opposite points on the top of the cone. We
                // set the cone angle to be that for lone pairs in an sp3 setting
                // and hi-jack the LonePairAddition code to generate two points.
                assert(lonePairs.size() == 1);
                const auto diff = coord - lonePairs[0];
                const auto backward = coord - diff;
                // in sp3 system angle between lp is 109.47
                double angle = 109.47 * M_PI / 180.0;
                // This is the size of the cone edge- so that cone length is the
                // same // as the lone pair distance.
                double dist = diff.length();
                RDGeom::Point3D point2, point3;
                AcceptorAtom::addTwoPairsRandomlyToTrigonal(coord, backward, point2,
                                                            point3, angle / 2, dist);
                return std::make_unique<ConeFeatureGeometry>(coord, point2, point3);
            }
            case HydrogenBondGeometry::None:
                assert(lonePairs.empty());
                return std::make_unique<SphereFeatureGeometry>(coord, hBondLen);
            default:
                throw std::domain_error("Unknown acceptor geometry");
        }
    }
} // Gape
