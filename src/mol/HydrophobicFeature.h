#pragma once

#include "Feature.h"

using namespace RDKit;

namespace Gape
{
	/**
	 * Represents a Hydrophobic atom (Carbon). This feature is used in the fitting
	 * process, but is not used in scoring.
	 * 
	 * @author Gareth Jones
	 * 
	 */
	class HydrophobicFeature : public Feature
	{
	public:
		explicit HydrophobicFeature(const int featureSetNumber) : Feature(FeatureType::HydrophobicAtom, featureSetNumber,
		                                                         "HYDROPHOBIC_ATOM", true)
		{
		}

		HydrophobicFeature(const int featureSetNumber, const SuperpositionMolecule* spMol,
		                   const Atom* featureAtom);

		static std::vector<std::shared_ptr<Feature>> findHydrophobicFeatures(
			const SuperpositionMolecule* superpositionMolecule);

		FeatureScore score(const Feature& otherFeature, const SuperpositionCoordinates& coordinates,
		             const SuperpositionCoordinates& otherCoordinates)  const override;

		std::string pharmLabel() const override;

		std::string info() const override;

		void calculateCoordinates(SuperpositionCoordinates& superpositionCoordinates) const override;

		std::unique_ptr<PharmFeatureGeometry> getPharmFeatureGeometry(
			const SuperpositionCoordinates& superpositionCoordinates) const override;

		const RDGeom::Point3D& getFittingPoint(const SuperpositionCoordinates& superpositionCoordinates) const override;
	};
} // namespace Gape
