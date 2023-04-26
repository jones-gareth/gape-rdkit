//
// Created by Gareth Jones on 3/23/2023.
//

#include <boost/format.hpp>
#include <boost/algorithm/string.hpp>
#include "AcceptorAtomFeature.h"

#include "HydrogenBondingType.h"
#include "../util/Reporter.h"

namespace Gape
{
	void AcceptorAtomFeature::getSolvationPoint(const SuperpositionCoordinates& superpositionCoordinates,
	                                            RDGeom::Point3D& solvationPoint) const
	{
		const auto& lonePairs = superpositionCoordinates.getFeatureCoordinates(FeatureType::AcceptorAtomFeature, atom);
		if (lonePairs.size() == 1)
		{
			solvationPoint = lonePairs[0];
		}
		else if (lonePairs.size() == 2)
		{
			const auto midPoint = (lonePairs[0] + lonePairs[1]) / 2.0;
			const auto& coordinate = superpositionCoordinates.getConformer().getAtomPos(atom->getIdx());
			auto diff = midPoint - coordinate;
			diff.normalize();
			diff *= hBondLen;
			solvationPoint = coordinate + diff;
		}
		else
		{
			assert(false);
		}
	}

	thread_local double AcceptorAtomFeature::hBondLen = 2.9;
	thread_local double AcceptorAtomFeature::chargeFactor = 2.0;
	thread_local double AcceptorAtomFeature::matchFactor = 1.0;
	thread_local bool AcceptorAtomFeature::scaleLonePairs = true;

	std::vector<std::shared_ptr<Feature>>
	AcceptorAtomFeature::findAcceptorAtoms(const SuperpositionMolecule* superpositionMolecule)
	{
		auto& mol = superpositionMolecule->getMol();
		auto& acceptors = superpositionMolecule->getAcceptors();
		std::vector<std::shared_ptr<Feature>> features;

		int featureSetNumber = 0;
		for (auto [atom, t] : acceptors)
		{
			auto acceptorAtom = mol.getAtomWithIdx(atom->getIdx());
			auto feature = std::make_shared<AcceptorAtomFeature>(featureSetNumber++, superpositionMolecule,
			                                                     acceptorAtom);
			features.push_back(std::static_pointer_cast<Feature>(feature));
		}
		return features;
	}

	AcceptorAtomFeature::AcceptorAtomFeature(const int featureSetNum, const SuperpositionMolecule* spMol,
	                                         const Atom* featureAtom): AcceptorAtomFeature(featureSetNum)
	{
		molecule = spMol;
		atom = featureAtom;
		const auto& acceptors = molecule->getAcceptors();
		const auto acceptorIt = acceptors.find(atom);
		assert(acceptorIt != acceptors.end());
		hydrogenBondingType = acceptorIt->second.get();
		charged = false; //TODO add charge to hydrogen bonding types
		REPORT(Reporter::DEBUG) << "N Lone Pairs " << nLonePairs;
		acceptorAtom = std::make_unique<AcceptorAtom>(molecule, atom);
	}

	std::string AcceptorAtomFeature::pharmLabel() const
	{
		auto name = boost::to_upper_copy<std::string>(hydrogenBondingType->name);
		std::replace(name.begin(), name.end(), ' ', '_');
		std::string options;
		if (charged)
		{
			options += " [charged = yes]";
		}
		const auto label = (boost::format("ATOM_%d") % (atom->getIdx() + 1)).str();
		return featureSetName + " " + label + options;
	}

	std::string AcceptorAtomFeature::info() const
	{
		const auto format = boost::format("Acceptor [%d]") % (atom->getIdx() + 1);
		return format.str();
	}

	const RDGeom::Point3D& AcceptorAtomFeature::getFittingPoint(
		const SuperpositionCoordinates& superpositionCoordinates) const
	{
		return superpositionCoordinates.getConformer().getAtomPos(atom->getIdx());
	}

	void AcceptorAtomFeature::calculateCoordinates(SuperpositionCoordinates& superpositionCoordinates) const
	{
		std::vector<RDGeom::Point3D> lonePairs;
		const auto& conformer = superpositionCoordinates.getConformer();
		acceptorAtom->addLonePairs(conformer, lonePairs);
		superpositionCoordinates.addFeatureCoordinates(FeatureType::AcceptorAtomFeature, atom, lonePairs);
	}

	double AcceptorAtomFeature::score(const Feature& otherFeature, const SuperpositionCoordinates& coordinates,
	                                  const SuperpositionCoordinates& otherCoordinates)
	{
		const auto& other = dynamic_cast<const AcceptorAtomFeature&>(otherFeature);

		// Score for acceptor overlay
		const auto& coordinate = coordinates.getConformer().getAtomPos(atom->getIdx());
		const auto& otherCoordinate = otherCoordinates.getConformer().getAtomPos(other.atom->getIdx());
		const auto sqrDistance = (coordinate - otherCoordinate).lengthSq();
		auto accScore = Feature::score(sqrDistance);
		if (accScore > maximumGaussianScore)
		{
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
		if (accScore <= 0)
			return 0;
	}
} // Gape
