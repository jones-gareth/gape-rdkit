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

		thread_local static double maxDonorAngle;
		thread_local static double minDonorAngle;
		thread_local static bool scoreDonorAtoms;

		bool charged = false;

	public:
		DonorHydrogenFeature(const int featureSetNumber) : Feature(FeatureType::DonorInteractionPoint, featureSetNumber,
		                                                           "DONOR_HYDROGEN", false)
		{
		};

		DonorHydrogenFeature(const int featureSetNumber, SuperpositionMolecule* spMol,
		                     Atom* atom): DonorHydrogenFeature(featureSetNumber)
		{
			molecule = spMol;
			this->atom = atom;
			const auto &mol = molecule->getMol();
			for (const auto nbr: mol.atomNeighbors(atom))
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

		static std::vector<const Atom*> findDonorHydrogens(const SuperpositionMolecule* molecule);
	};
}
