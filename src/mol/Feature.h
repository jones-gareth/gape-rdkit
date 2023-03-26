//
// Created by Gareth Jones on 1/19/23.
//

#ifndef GAPE_FEATURE_H
#define GAPE_FEATURE_H

#include <GraphMol/GraphMol.h>

#include "../gape/SuperpositionMolecule.h"
#include "PharmFeatureGeometry.h"

using namespace RDKit;

namespace Gape
{
	enum FeatureType
	{
		HydrophobicAtom,
		DonorInteractionPoint,
		AcceptorAtomFeature,
		AromaticRing
		// TODO UserFeatures
	};

	class Feature
	{
	protected:
		const FeatureType featureType;
		/**
		 * This is the feature number- e.g. HYDROPHOBIC_ATOM number 2
		 */
		const int featureSetNumber;

		/**
		 * A descriptive name for the feature
		 */
		const std::string featureSetName;

		/**
		 * Set this if the feature is atom-centered, e.g. an acceptor
		 */
		const bool atomFeature;

		/**
		 * Atom if feature is atom centered
		 */
		const Atom* atom = nullptr;

		/**
		 * Molecule containing this feature
		 */
		const SuperpositionMolecule* molecule;

		/**
		 * Set this if we've stored a mapping
		 */
		bool mapped;

		/**
		 * If this feature (in the base molecule) is mapped to a feature in another
		 * molecule, we can store it here.
		 */
		Feature* mappedFeature;

		// sqr distance to other feature
		double squaredDistance;


		// These variables are used to track the best mappings so far- for
		// resolving one-to-one mappings.
		double bestScore, bestGeometricScore;

		/**
		 * If this feature (in another molecule) is mappped to a feature in the base
		 * molecule, we can store the base molecule feature here.
		 */
		Feature* baseFeature;

		/**
		 * A feature normally has an point which can be considered it's interaction
		 * point. E.g. for acceptors the atom, for donors the virtual acceptor on
		 * the donor-hydrogen line and for rings the ring center. We can use this
		 * point in chromosome encoding for mapping features and place Gaussians on
		 * this point for scoring feature overlap.
		 */
		RDGeom::Point3D coordinate;

		/**
		 * Set this true if we consider this feature to be part of a pharmacophore.
		 */
		bool pharmacophorePoint;

		/**
		 * number of molecules significantly contribution to this feature if it is a
		 * pharmacophore point.
		 *
		 */
		int numberMatched;

		/**
		 * Set this true if we can use this feature in chromosome encoding.
		 */
		bool mappingFeature = true;

		thread_local static double maximumGaussianScore;
		thread_local static double solventVolOk;

		/**
		 * The score function must set this variable
		 */
		double geometricScore;

		Feature(const FeatureType featureType, const int featureSetNumber, const std::string& featureSetName,
		        const bool atomFeature):
			featureType(featureType), featureSetNumber(featureSetNumber), featureSetName(featureSetName),
			atomFeature(atomFeature)
		{
			molecule = nullptr;
			mapped = false;
			mappedFeature = nullptr;
			squaredDistance = DBL_MAX;
			bestGeometricScore = -DBL_MAX;
			bestScore = -DBL_MAX;
			geometricScore = -DBL_MAX;
			pharmacophorePoint = false;
			baseFeature = nullptr;
			numberMatched = 0;
		};

		/**
		 * The alpha used in a gaussian to determine if a point is solvent
		 * accessible.
		 */
		static double solvationAlpha;

	private:
		// TODO move scoring and book-keeping to another class
		/**
		 * Gaussian parameters for the core built in features.
		 */
		thread_local static double alpha, gaussianN, radius;

	public:
		static const int N_BUILTIN_FEATURES = 4, MAX_FEATURES = 14;

		static void initialize();

		/**
		 * Geometry of feature
		 */
		std::unique_ptr<PharmFeatureGeometry> pharmFeatureGeometry;

		/**
		 * Each subclass must implement the score method which returns the
		 * similarity between two features.
		 *
		 * @param f
		 * @return
		 */
		virtual double score(const Feature& otherFeature, const Conformer& conformer, const Conformer& otherConformer);

		/**
		 * Returns a label for the feature to use in mol2 or sdf pharmacophore
		 * files.
		 *
		 * @return
		 */
		std::string atomLabel() const;

		/**
		 * Generates a label to use to describe the pharmacophore and feature
		 * geometry in the descriptions that appear in structure files (either in
		 * sd-field or a comment section in mol2 files). Similar to pharmLabel, but
		 * includes feature geometry definition
		 *
		 * @return
		 *
		 * @see #pharmLabel()
		 */
		std::string featureLabel() const;

		/**
		 * Generates a label to use to describe the pharmacophore in the
		 * descriptions that appear in structure files (either in sd-field or a
		 * comment section in mol2 files).
		 *
		 * @return
		 */
		virtual std::string pharmLabel() const;


		/**
		 * Returns an information/description string for this pharmacophore.
		 *
		 * @return
		 */
		virtual std::string info() const;

		/**
		 * Sets the gaussian radius to r. Adjusts the other gausian parameters so
		 * that the overlap is normalized. This is used only by the built-in
		 * features and not the user-defined features, which have their own defined
		 * radius.
		 *
		 * @param r
		 */
		static void setRadius(double r);

		/**
		 * Determines the integral between two gaussians. The square distance
		 * between gaussian centers is given. The gaussian parameter for alpha is
		 * taken for the class variable, so this method only works for the class
		 * features.
		 *
		 * @param sqrDistance
		 * @return
		 */
		static double score(double sqrDistance);


		/**
		 * Determines the integral between two gaussians. The square distance
		 * between gaussian centers is given as is the alpha parameter for the
		 * gaussians (both gaussians have the same alpha parameter).
		 *
		 * @param sqrDistance
		 * @param a
		 * @return
		 */
		static double score(double sqrDistance, double a);


		/**
		 * Determine a gaussian alpha parameter from a "radius". It's more intuitive
		 * to specify a radius rather that alpha.
		 *
		 * @param r
		 * @return
		 */
		static double getAlpha(double r);


		/**
		 * Set the alpha used in a gaussian to determine if a point is solvent
		 * accessible.
		 *
		 * @param r
		 */
		static void setSolvationAlpha(double r);

		/**
		 * Returns a penalty if point is not solvent accessible. Used to determine
		 * if acceptor lone pairs or donor hydrogen fitting points are available. A
		 * gaussian of alpha solvationAlpha is placed at point and we check for
		 * overlap with any atoms in mol. Overlap with a specific atom is ignored
		 * (this atom being the acceptor or donor hydrogen).
		 *
		 * @param point
		 * @param mol
		 * @param atom
		 * @return
		 */
		double solvationPenalty(const RDGeom::Point3D& point, const ROMol& mol, const Conformer& conformer,
		                        const Atom& atom) const;


		/**
		 * Returns the atom-centered coordinate. As a side effect copies the
		 * coordinate from the molecule to this class- if we do this other methods
		 * are free to manipulate the molecule coordinates without affecting the
		 * feature coordinate. For features without atom-centered coordinates
		 * (rings, donor hydrogens and some user-features) you need to override this
		 * method.
		 *
		 * @return
		 */
		virtual const RDGeom::Point3D& calculateCoordinate(const Conformer& conformer);

		virtual const PharmFeatureGeometry& getPharmFeatureGeometry(const Conformer& conformer);

		/**
		 * As getCoordinate, but returns the saved class copy of the coordinate.
		 *
		 * @return
		 */
		const RDGeom::Point3D& getSavedCoordinate() const;

		/**
		 * Returns the square distance between this feature and the other feature
		 * mapped to it.
		 */
		double calculateSqrDist();

		std::string mappingInfo() const;
	};
} // Gape

#endif //GAPE_FEATURE_H
