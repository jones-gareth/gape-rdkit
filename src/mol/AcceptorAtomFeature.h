//
// Created by Gareth Jones on 3/23/2023.
//

#pragma once

#include "Feature.h"
#include "HydrogenBondingType.h"
#include "AcceptorAtom.h"

using namespace RDKit;

namespace Gape
{
	class AcceptorAtomFeature : public Feature
	{
	private:
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
		/**
		 * For scoring functions based on lone pair overlap we can scale totals
		 * based on the number of lone-pairs
		 */
		thread_local static bool scaleLonePairs;

		std::vector<RDGeom::Point3D> lonePairCoords;

		int nLonePairs = 0;

		// Coordinate for storing a point that should be solvent accessible
		RDGeom::Point3D solvationPoint;

		// Normal coordinate- for planar acceptors
		RDGeom::Point3D normal;


		// The default is to use Mills and Dean types in preference to GASP/GOLD
		// I've now stopped using the GASP/GOLD types completely

		// Acceptor geometry e.g. LP, PLANE
		HydrogenBondGeometry geometry = HydrogenBondGeometry::None;

		/**
		 * Set for a charged acceptor
		 */
		bool charged = false;

		std::unique_ptr<Gape::AcceptorAtom> acceptorAtom;

	public:
		AcceptorAtomFeature(const int featureSetNumber) : Feature(FeatureType::AcceptorAtomFeature, featureSetNumber,
		                                                          "ACCEPTOR_ATOM", true)
		{
		}


		AcceptorAtomFeature(const int featureSetNum, const SuperpositionMolecule* spMol,
		                    const Atom* featureAtom);

		static std::vector<std::shared_ptr<Feature>> findAcceptorAtoms(
			const SuperpositionMolecule* superpositionMolecule);
	};
} // Gape
