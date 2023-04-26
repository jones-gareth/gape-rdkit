//
// Created by Gareth Jones on 3/23/2023.
//

#pragma once

#include "Feature.h"
#include "AcceptorAtom.h"

using namespace RDKit;

namespace Gape
{
	class HydrogenBondingType;
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

		// std::vector<RDGeom::Point3D> lonePairCoords;

		int nLonePairs = 0;

		// Coordinate for storing a point that should be solvent accessible
		RDGeom::Point3D solvationPoint;

		// Normal coordinate- for planar acceptors
		RDGeom::Point3D normal;


		// The default is to use Mills and Dean types in preference to GASP/GOLD
		// I've now stopped using the GASP/GOLD types completely

		// Acceptor geometry e.g. LP, PLANE
		// HydrogenBondGeometry geometry = HydrogenBondGeometry::None;

		/**
		 * Set for a charged acceptor
		 */
		bool charged = false;

		std::unique_ptr<Gape::AcceptorAtom> acceptorAtom;

		/**
		 * Determines a point that should be solvent accessible if an acceptor with
		 * one or two lone pairs is able to accept a hydrogen bond. Note that in the
		 * acceptor geometry model an acceptor with three sp3 lone pairs has a
		 * single lone pair in the forward direction.
		 */
		void getSolvationPoint(const SuperpositionCoordinates& superpositionCoordinates,
		                       RDGeom::Point3D& solvationPoint) const;


	public:
		AcceptorAtomFeature(const int featureSetNumber) : Feature(FeatureType::AcceptorAtomFeature, featureSetNumber,
		                                                          "ACCEPTOR_ATOM", true)
		{
		}


		AcceptorAtomFeature(const int featureSetNum, const SuperpositionMolecule* spMol,
		                    const Atom* featureAtom);

		static std::vector<std::shared_ptr<Feature>> findAcceptorAtoms(
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
} // Gape
