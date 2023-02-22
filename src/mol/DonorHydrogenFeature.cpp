#include "DonorHydrogenFeature.h"

#include <boost/format.hpp>
#include <boost/algorithm/string.hpp>

#include "HydrogenBondingType.h"

namespace Gape
{
	thread_local double DonorHydrogenFeature::hBondLen = 2.9;
	thread_local double DonorHydrogenFeature::chargeFactor = 2.0;

	thread_local double DonorHydrogenFeature::maxDonorAngle = 120.0 * M_PI / 180.0;
	thread_local double DonorHydrogenFeature::minDonorAngle = 90 * M_PI / 180.0;
	thread_local bool DonorHydrogenFeature::scoreDonorAtoms = true;

	std::vector<std::shared_ptr<Feature>> DonorHydrogenFeature::findDonorHydrogens(
		const SuperpositionMolecule* superpositionMolecule)
	{
		auto& mol = superpositionMolecule->getMol();
		auto& donors = superpositionMolecule->getDonors();
		std::vector<std::shared_ptr<Feature>> features;
		int featureSetNumber = 0;
		for (const auto atom : mol.atoms())
		{
			if (const auto donorIt = donors.find(atom); donorIt != donors.end())
			{
				for (const auto neighbor : mol.atomNeighbors(atom))
				{
					if (neighbor->getAtomicNum() == 1)
					{
						auto feature = std::make_shared<DonorHydrogenFeature>(
							featureSetNumber++, superpositionMolecule, neighbor);
						features.push_back(std::static_pointer_cast<Feature>(feature));
					}
				}
			}
		}
		return features;
	}

	DonorHydrogenFeature::DonorHydrogenFeature(const int featureSetNum, const SuperpositionMolecule* spMol,
	                                           Atom* featureAtom) : DonorHydrogenFeature(featureSetNum)
	{
		molecule = spMol;
		atom = featureAtom;
		const auto& mol = molecule->getMol();
		for (const auto nbr : mol.atomNeighbors(atom))
		{
			assert(donor == nullptr);
			donor = nbr;
		}
		const auto& donors = molecule->getDonors();
		const auto donorIt = donors.find(donor);
		assert(donorIt != donors.end());
		hydrogenBondingType = donorIt->second.get();
		charged = false; //TODO add charge to hydrogen bonding types
	}

	RDGeom::Point3D& DonorHydrogenFeature::calculateCoordinate(const Conformer& conformer)
	{
		const auto& donorCoord = conformer.getAtomPos(donor->getIdx());
		const auto& hydrogenCoord = conformer.getAtomPos(atom->getIdx());
		auto bondVector = donorCoord.directionVector(hydrogenCoord);
		bondVector *= hBondLen;
		coordinate = donorCoord + bondVector;
		return coordinate;
	}

	std::string DonorHydrogenFeature::info() const
	{
		const int atomNum = atom->getIdx() + 1;
		const auto format = boost::format("Donor H [%s]") % atomNum;
		return format.str();
	}

	std::string DonorHydrogenFeature::pharmLabel() const
	{
		auto name = boost::to_upper_copy<std::string>(hydrogenBondingType->name);
		std::replace(name.begin(), name.end(), ' ', '_');
		std::string options("");
		if (charged)
		{
			options += " [charged = yes]";
		}
		const auto label = (boost::format("ATOM_%d") % (atom->getIdx() + 1)).str();
		return featureSetName + " " + label + options;
	}

	double DonorHydrogenFeature::score(const Feature &f)
{
			const auto &other = dynamic_cast<const DonorHydrogenFeature&>(f);

			// fitting point Gaussian
			coordinate.
			double vol = Feature.score(Coord.sqrDistance(coordinate, other.coordinate));
			Coord.midPoint(coordinate, other.coordinate, midPoint.get());

			double correct = .0;
			if (!virtual) {
				// Corrections to make sure fitting point is solvent accessible
				double molVol = Feature.solvationPenalty(midPoint.get(), molecule, atom);
				double otherMolVol = Feature.solvationPenalty(midPoint.get(), other.molecule,
						other.atom);

				if (molVol < .0)
					molVol = .0;
				if (otherMolVol < .0)
					otherMolVol = .0;
				correct = molVol > otherMolVol ? molVol : otherMolVol;

				if (logDebug) {
					logger.debug(info() + " " + other.info() + " fitting vol " + vol
							+ " molA " + molVol + " molB " + otherMolVol);
				}
			}
			double fitScore = vol - correct;
			if (fitScore < 0)
				fitScore = 0;

			geometricScore = fitScore;

			double score;
			if (scoreDonorAtoms.get()) {
				// now check that donor atoms are reasonably positioned.
				if (fitScore == 0)
					return .0;
				Coord.subtract(donorCoord, midPoint.get(), vec1.get());
				Coord.subtract(other.donorCoord, midPoint.get(), vec2.get());
				double angle = Coord.angle(vec1.get(), vec2.get());
				if (logDebug)
					logger.debug("Angle " + angle + " max " + maxDonorDonorAngle + " min "
							+ minDonorDonorAngle);
				double donorScore = .0;
				if (angle > maxDonorDonorAngle.get())
					return .0;
				else if (angle < minDonorDonorAngle.get())
					donorScore = 1.0;
				else {
					donorScore = 1 - (angle - minDonorDonorAngle.get())
						/ (maxDonorDonorAngle.get() - minDonorDonorAngle.get());
				}

				if (logDebug)
					logger.debug(" donor score " + donorScore);
				score = fitScore * donorScore;
			}
			else {
				score = fitScore;
				if (fitScore > maximumGaussianScore.get())
					score = maximumGaussianScore.get();
			}

			if (logDebug)
				logger.debug(" geometric score " + geometricScore);

			if (score < .0)
				return .0;

			// Mills and Dean probabilities for types
			double prob = hydrogenBondingType.getProbability()
				+ other.hydrogenBondingType.getProbability();
			score *= prob;

			// increase for matching types
			if (hydrogenBondingType.getId() == other.hydrogenBondingType.getId())
				score *= matchFactor.get();

			// increase if both donors are charged
			if (charged && other.charged)
				score *= chargeFactor.get();

			if (logDebug)
				logger.debug(" type score " + score);
			return score;
		}

}
