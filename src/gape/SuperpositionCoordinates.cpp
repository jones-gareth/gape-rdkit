#include "SuperpositionCoordinates.h"

namespace Gape
{
	void SuperpositionCoordinates::addFeatureCoordinates(FeatureType featureType, const Atom* atom,
	                                                     const std::vector<RDGeom::Point3D>& coordinates)
	{
		if (featureCoordinates.find(featureType) == featureCoordinates.end())
		{
			const std::map<const Atom*, std::vector<RDGeom::Point3D>> featureTypeCoordinates;
			featureCoordinates.emplace(featureType, featureTypeCoordinates);
		}
		auto featureTypeCoordinates = featureCoordinates[featureType];
		featureTypeCoordinates.emplace(atom, coordinates);
	}
}
