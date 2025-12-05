#include "SuperpositionCoordinates.h"
#include <Geometry/Transform3D.h>
#include "SuperpositionMolecule.h"

namespace Gape
{
	void SuperpositionCoordinates::addFeatureCoordinates(FeatureType featureType, const Atom* atom,
	                                                     const std::vector<RDGeom::Point3D>& coordinates)
	{
		assert(&atom->getOwningMol() == &molecule->getMol());
		if (featureCoordinates.find(featureType) == featureCoordinates.end())
		{
			const std::map<const Atom*, std::vector<RDGeom::Point3D>> featureTypeCoordinates;
			featureCoordinates.emplace(featureType, featureTypeCoordinates);
		}
		auto  &featureTypeCoordinates = featureCoordinates[featureType];
		featureTypeCoordinates.emplace(atom, coordinates);
	}

	const std::vector<RDGeom::Point3D> &SuperpositionCoordinates::getFeatureCoordinates(FeatureType featureType, const Atom *atom) const {
		assert(&atom->getOwningMol() == &molecule->getMol());
		return featureCoordinates.at(featureType).at(atom);
	}

	std::vector<const Atom*> SuperpositionCoordinates::getFeatureAtoms(const FeatureType featureType) const
	{
		std::vector<const Atom*> atoms;
		if (featureCoordinates.find(featureType) == featureCoordinates.end())
		{
			return atoms;
		}

		const auto& featureTypeCoordinates = featureCoordinates.at(featureType);
		for (const auto& [fst, snd]: featureTypeCoordinates)
		{
			atoms.push_back(fst);
		}
		return atoms;
	}

    void SuperpositionCoordinates::transformCoordinates(const Eigen::Transform<double, 3, Eigen::Affine> &transform) {
        RDGeom::Transform3D rdTransform;
        for (auto row=0; row<4; ++row) {
            for (auto col=0; col<4; col++) {
                rdTransform.setVal(row, col, transform(row, col));
            }
        }
        for (auto &point: conformer.getPositions()) {
            rdTransform.TransformPoint(point);
        }
        for (auto &featureTypeCoordinatesEntry: featureCoordinates) {
            for (auto &atomEntry: featureTypeCoordinatesEntry.second) {
                for (auto &featureCoordinate: atomEntry.second) {
                    rdTransform.TransformPoint(featureCoordinate);
                }
            }
        }
    }

    /*
    const RDGeom::Point3D &SuperpositionCoordinates::getFittingCoordinate(const Gape::Feature &feature) const {
        if (feature.atomFeature) {
            const auto index = feature.getAtom()->getIdx();
            return conformer.getAtomPos(index);
        }

        const auto& featureTypeCoordinates = featureCoordinates.at(feature.getFeatureType());
        const auto& coords = featureTypeCoordinates.at(feature.getAtom());
        return coords[0];
    }
     */
}
