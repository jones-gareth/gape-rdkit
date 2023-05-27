#pragma once

#include "Feature.h"

using namespace RDKit;

namespace Gape
{
	class AromaticRingFeature : public Feature
	{
		std::vector<const Atom*> ringAtoms;

		// Strict aromatic ring
		bool aromatic;

		// Planar Ring
		bool planar;


		thread_local static double normalLength;

	public:
		AromaticRingFeature(int featureSetNumber) :
			Feature(FeatureType::AromaticRing, featureSetNumber, "Aromatic Ring", false)
		{
		}

		AromaticRingFeature(const int featureSetNumber, const SuperpositionMolecule* spMol,
		                    const std::vector<const Atom*>& ringAtoms, const bool aromatic, const bool planar):
			AromaticRingFeature(featureSetNumber)
		{
			this->ringAtoms = ringAtoms;
			this->aromatic = aromatic;
			this->planar = planar;
			molecule = spMol;
		}

		static std::vector<std::shared_ptr<Feature>> findAromaticRings(
			const SuperpositionMolecule* superpositionMolecule);

		double score(const Feature& otherFeature, const SuperpositionCoordinates& coordinates,
		             const SuperpositionCoordinates& otherCoordinates) override;
		std::string pharmLabel() const override;
		std::string info() const override;
		void calculateCoordinates(SuperpositionCoordinates& superpositionCoordinates) const override;
		std::unique_ptr<PharmFeatureGeometry> getPharmFeatureGeometry(
			const SuperpositionCoordinates& superpositionCoordinates) const override;
		const RDGeom::Point3D& getFittingPoint(const SuperpositionCoordinates& superpositionCoordinates) const override;
	};
} // namespace Gape
