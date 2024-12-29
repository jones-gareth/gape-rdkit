#include "HydrophobicFeature.h"
#include "util/Reporter.h"
#include "gape/SuperpositionMolecule.h"
#include <boost/format.hpp>
#include <GraphMol/GraphMol.h>
#include <GraphMol/PeriodicTable.h>

namespace Gape
{
	HydrophobicFeature::HydrophobicFeature(const int featureSetNumber, const SuperpositionMolecule* spMol,
	                                       const Atom* featureAtom) : Feature(
		FeatureType::DonorInteractionPoint, featureSetNumber,
		"HYDROPHOBIC_ATOM", true)
	{
	}


	std::vector<std::shared_ptr<Feature>> HydrophobicFeature::findHydrophobicFeatures(
		const SuperpositionMolecule* superpositionMolecule)
	{
		auto& mol = superpositionMolecule->getMol();
		std::vector<std::shared_ptr<Feature>> features;

		int featureSetNumber = 0;
		for (const auto atom : mol.atoms())
		{
			{
				if (atom->getAtomicNum() == 6 && atom->getHybridization() != Atom::SP)
				{
					auto acceptorAtom = mol.getAtomWithIdx(atom->getIdx());
					auto feature = std::make_shared<HydrophobicFeature>(featureSetNumber++, superpositionMolecule,
					                                                    acceptorAtom);
					features.push_back(std::static_pointer_cast<Feature>(feature));
				}
			}
		}
		return features;
	}

	FeatureScore HydrophobicFeature::score(const Feature& otherFeature, const SuperpositionCoordinates& coordinates,
	                                 const SuperpositionCoordinates& otherCoordinates) const
	{
		const FeatureScore result(.0, .0);
		return result;
	}

	std::string HydrophobicFeature::pharmLabel() const
	{
		const auto label = (boost::format("ATOM_%d") % (atom->getIdx() + 1)).str();
		return featureSetName + " " + label;
	}

	std::string HydrophobicFeature::info() const
	{
		const auto format = boost::format("HydroP Atom [%d]") % (atom->getIdx() + 1);
		return format.str();
	}

	void HydrophobicFeature::calculateCoordinates(SuperpositionCoordinates& superpositionCoordinates) const
	{
	}

	std::unique_ptr<PharmFeatureGeometry> HydrophobicFeature::getPharmFeatureGeometry(
		const SuperpositionCoordinates& superpositionCoordinates) const
	{
		const auto radius = PeriodicTable::getTable()->getRvdw(6);
		const auto& coord = superpositionCoordinates.getConformer().getAtomPos(atom->getIdx());
		return std::make_unique<SphereFeatureGeometry>(coord, radius);
	}

	const RDGeom::Point3D& HydrophobicFeature::getFittingPoint(
		const SuperpositionCoordinates& superpositionCoordinates) const
	{
		return superpositionCoordinates.getConformer().getAtomPos(atom->getIdx());
	}
} //namespace Gape
