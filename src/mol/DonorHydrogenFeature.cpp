#include "DonorHydrogenFeature.h"

#include <boost/format.hpp>
#include <boost/algorithm/string.hpp>

#include "HydrogenBondingType.h"
#include "util/Reporter.h"
#include "gape/SuperpositionMolecule.h"
#include "mol/PartialCharge.h"

namespace Gape {
    thread_local double DonorHydrogenFeature::hBondLen = 2.9;
    thread_local double DonorHydrogenFeature::chargeFactor = 2.0;
    thread_local double DonorHydrogenFeature::matchFactor = 1.0;

    thread_local double DonorHydrogenFeature::maxDonorDonorAngle = 100.0 * M_PI / 180.0;
    thread_local double DonorHydrogenFeature::minDonorDonorAngle = 60 * M_PI / 180.0;
    thread_local bool DonorHydrogenFeature::scoreDonorAtoms = true;

    double DonorHydrogenFeature::minDonorPartialCharge = 0.2;

    std::vector<std::shared_ptr<Feature> > DonorHydrogenFeature::findDonorHydrogens(
        const SuperpositionMolecule *superpositionMolecule) {
        auto &mol = superpositionMolecule->getMol();
        auto &donors = superpositionMolecule->getDonors();
        std::vector<std::shared_ptr<Feature> > features;
        int featureSetNumber = 0;
        for (const auto atom: mol.atoms()) {
            if (const auto donorIt = donors.find(atom); donorIt != donors.end()) {
                for (const auto neighbor: mol.atomNeighbors(atom)) {
                    if (neighbor->getAtomicNum() == 1) {
                        auto feature = std::make_shared<DonorHydrogenFeature>(
                            featureSetNumber++, superpositionMolecule, neighbor);
                        features.push_back(std::static_pointer_cast<Feature>(feature));
                    }
                }
            }
        }
        return features;
    }

    DonorHydrogenFeature::DonorHydrogenFeature(const int featureSetNum, const SuperpositionMolecule *spMol,
                                               const Atom *featureAtom) : DonorHydrogenFeature(featureSetNum) {
        molecule = spMol;
        atom = featureAtom;
        const auto &mol = molecule->getMol();
        for (const auto nbr: mol.atomNeighbors(atom)) {
            assert(donor == nullptr);
            donor = nbr;
        }
        const auto &donors = molecule->getDonors();
        const auto donorIt = donors.find(donor);
        assert(donorIt != donors.end());
        hydrogenBondingType = donorIt->second.get();
        charged = PartialCharge::getPartialCharge(donor) > minDonorPartialCharge;
    }

    void DonorHydrogenFeature::calculateCoordinates(SuperpositionCoordinates &superpositionCoordinates) const {
        const auto &conformer = superpositionCoordinates.getConformer();
        const auto &donorCoord = conformer.getAtomPos(donor->getIdx());
        const auto &hydrogenCoord = conformer.getAtomPos(atom->getIdx());
        auto bondVector = donorCoord.directionVector(hydrogenCoord);
        bondVector *= hBondLen;
        const std::vector<RDGeom::Point3D> coordinates{donorCoord + bondVector};
        superpositionCoordinates.addFeatureCoordinates(FeatureType::DonorInteractionPoint, atom, coordinates);
    }

    std::string DonorHydrogenFeature::info() const {
        const int atomNum = atom->getIdx() + 1;
        const auto format = boost::format("Donor H [%s]") % atomNum;
        return format.str();
    }

    std::string DonorHydrogenFeature::pharmLabel() const {
        auto name = boost::to_upper_copy<std::string>(hydrogenBondingType->name);
        std::replace(name.begin(), name.end(), ' ', '_');
        std::string options;
        if (charged) {
            options += " [charged = yes]";
        }
        const auto label = (boost::format("ATOM_%d") % (atom->getIdx() + 1)).str();
        return featureSetName + " " + label + options;
    }

    FeatureScore DonorHydrogenFeature::score(const Feature &otherFeature, const SuperpositionCoordinates &coordinates,
                                             const SuperpositionCoordinates &otherCoordinates) const {
        const auto &other = dynamic_cast<const DonorHydrogenFeature &>(otherFeature);

        // fitting point Gaussian
        const auto &coordinate = coordinates.getFeatureCoordinates(FeatureType::DonorInteractionPoint, atom)[0];
        const auto &otherCoordinate = otherCoordinates.getFeatureCoordinates(
            FeatureType::DonorInteractionPoint, other.atom)[0];
        const auto sqrDistance = (coordinate - otherCoordinate).lengthSq();
        const auto vol = Feature::score(sqrDistance);
        // TODO - check for maximum Gaussian score?
        const auto midPoint = (coordinate + otherCoordinate) / 2.0;
        const auto mol = molecule->getMol();
        const auto otherMol = other.molecule->getMol();

        // Corrections to make sure fitting point is solvent accessible
        auto molVol = solvationPenalty(midPoint, coordinates);
        auto otherMolVol = otherFeature.solvationPenalty(midPoint, otherCoordinates);

        if (molVol < .0)
            molVol = .0;
        if (otherMolVol < .0)
            otherMolVol = .0;
        const auto correct = molVol > otherMolVol ? molVol : otherMolVol;

        REPORT(Reporter::DEBUG) <<
				info() << " " << other.info() << " fitting vol " << vol
					<< " molA " << molVol << " molB " << otherMolVol;

        double fitScore = vol - correct;
        if (fitScore < 0)
            fitScore = 0;

        const auto geometricScore = fitScore;
        double score;
        if (scoreDonorAtoms) {
            // now check that donor atoms are reasonably positioned.
            if (fitScore == 0) {
                FeatureScore result(.0, .0);
                return result;
            }
            const auto vec1 = donorCoordinate - midPoint;
            const auto vec2 = other.donorCoordinate - midPoint;
            const double angle = vec1.angleTo(vec2);
            REPORT(Reporter::DEBUG) << "Angle " << angle << " max " << maxDonorDonorAngle << " min "
					<< minDonorDonorAngle;
            double donorScore = .0;
            if (angle > maxDonorDonorAngle) {
                FeatureScore result(.0, .0);
                return result;
            }
            if (angle < minDonorDonorAngle) {
                donorScore = 1.0;
            } else {
                donorScore = 1 - (angle - minDonorDonorAngle)
                             / (maxDonorDonorAngle - minDonorDonorAngle);
            }
            REPORT(Reporter::DEBUG) << " donor score " << donorScore;
            score = fitScore * donorScore;
        } else {
            score = fitScore;
            if (fitScore > maximumGaussianScore) {
                score = maximumGaussianScore;
            }
        }

        REPORT(Reporter::DEBUG) << " geometric score " << geometricScore;

        if (score < .0) {
            FeatureScore result(.0, .0);
            return result;
        }

        // Mills and Dean probabilities for types
        double prob = hydrogenBondingType->probability
                      + other.hydrogenBondingType->probability;
        score *= prob;

        // increase for matching types
        if (hydrogenBondingType->name == other.hydrogenBondingType->name) {
            score *= matchFactor;
        }

        // increase if both donors are charged
        if (charged && other.charged) {
            score *= chargeFactor;
        }

        REPORT(Reporter::DEBUG) << " type score " << score;
        const FeatureScore result(score, geometricScore);
        return result;
    }

    std::unique_ptr<PharmFeatureGeometry> DonorHydrogenFeature::getPharmFeatureGeometry(
        const SuperpositionCoordinates &superpositionCoordinates) const {
        const auto &conformer = superpositionCoordinates.getConformer();
        auto &point1 = conformer.getAtomPos(donor->getIdx());
        auto &point2 = superpositionCoordinates.getFeatureCoordinates(FeatureType::DonorInteractionPoint, atom)[0];
        auto pharmFeatureGeometry = std::make_unique<VectorPharmFeatureGeometry>(point1, point2);
        return pharmFeatureGeometry;
    }

    const RDGeom::Point3D &DonorHydrogenFeature::getFittingPoint(
        const SuperpositionCoordinates &superpositionCoordinates) const {
        return superpositionCoordinates.getFeatureCoordinates(FeatureType::DonorInteractionPoint, atom)[0];
    }
}
