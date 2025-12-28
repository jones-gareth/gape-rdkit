#pragma once

#include "Feature.h"

using namespace RDKit;

namespace Gape
{
	class AromaticRingFeature : public Feature
	{
		std::vector<const Atom*> ringAtoms;

		// Strict aromatic ring
		bool aromatic = false;

		// Planar Ring
		bool planar = false;

		const static double normalLength;

	public:
		AromaticRingFeature(int featureSetNumber) :
			Feature(FeatureType::AromaticRing, featureSetNumber, "Aromatic Ring", false)
		{
		}

		AromaticRingFeature(int featureSetNumber, const SuperpositionMolecule* spMol,
							const std::vector<const Atom*>& ringAtoms, bool aromatic, bool planar, const std::vector<const Atom *> &atomsInUse);


		static std::vector<std::shared_ptr<Feature>> findAromaticRings(
			const SuperpositionMolecule* superpositionMolecule);

		FeatureScore score(const Feature& otherFeature, const SuperpositionCoordinates& coordinates,
		             const SuperpositionCoordinates& otherCoordinates) const override;
		std::string pharmLabel() const override;
		std::string info() const override;
		void calculateCoordinates(SuperpositionCoordinates& superpositionCoordinates) const override;
		std::unique_ptr<PharmFeatureGeometry> getPharmFeatureGeometry(
			const SuperpositionCoordinates& superpositionCoordinates) const override;
		const RDGeom::Point3D& getFittingPoint(const SuperpositionCoordinates& superpositionCoordinates) const override;

		RDGeom::Point3D ringCenter(const SuperpositionCoordinates& superpositionCoordinates) const;
		std::vector<RDGeom::Point3D> ringNormals(const SuperpositionCoordinates& superpositionCoordinates) const;
	};
} // namespace Gape
