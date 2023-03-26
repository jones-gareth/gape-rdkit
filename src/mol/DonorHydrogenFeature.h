#pragma once

#include "Feature.h"

using namespace RDKit;

namespace Gape
{
	class DonorHydrogenFeature : public Feature
	{
	private:
		const Atom* donor = nullptr;
		RDGeom::Point3D donorCoordinate;
		const HydrogenBondingType* hydrogenBondingType = nullptr;

		/**
		 * Length of a hydrogen bond. May be changed.
		 */
		thread_local static double hBondLen;

		/**
		 * Pairwise feature score is scaled by this if both features have negative
		 * charge. May be changed.
		 */
		thread_local static double chargeFactor;
		/**
		 * Pairwise feature score is scaled by this if both features have the same
		 * type. May be changed.
		 */
		thread_local static double matchFactor;

		thread_local static double maxDonorDonorAngle;
		thread_local static double minDonorDonorAngle;
		thread_local static bool scoreDonorAtoms;

		bool charged = false;

	public:
		DonorHydrogenFeature(const int featureSetNumber) : Feature(FeatureType::DonorInteractionPoint, featureSetNumber,
		                                                           "DONOR_HYDROGEN", false)
		{
		}

		DonorHydrogenFeature(const int featureSetNum, const SuperpositionMolecule* spMol,
		                     const Atom* featureAtom);

		static std::vector<std::shared_ptr<Feature>> findDonorHydrogens(
			const SuperpositionMolecule* superpositionMolecule);


		/*
		 * Determines the fitting point for the donor hydrogen. (non-Javadoc)
		 * 
		 * @see com.cairn.gape.Feature#getCoordinate()
		 */
		RDGeom::Point3D& calculateCoordinate(const Conformer& conformer) override;


		/*
		 * Descriptive string for the feature.
		 * 
		 * @see com.cairn.gape.Feature#info()
		 */
		std::string info() const override;

		std::string pharmLabel() const override;


		/*
		 * Returns similarity between two donor hydrogens. 
		 *
		 */
		double score(const Feature &otherFeature, const Conformer &conformer, const Conformer &otherConformer) override;


		const PharmFeatureGeometry& getPharmFeatureGeometry(const Conformer& conformer) override;
	};
}
