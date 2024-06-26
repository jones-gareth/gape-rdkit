﻿#pragma once

#include <GraphMol/GraphMol.h>
#include "mol/FeatureType.h"
#include <Eigen/Dense>

using namespace RDKit;

namespace Gape
{

	class SuperpositionCoordinates
	{
		Conformer conformer;
		std::map<FeatureType, std::map<const Atom*, std::vector<RDGeom::Point3D>>> featureCoordinates;

	public:
		explicit SuperpositionCoordinates(Conformer c) : conformer(std::move(c))
		{
		}

		void addFeatureCoordinates(FeatureType featureType, const Atom* atom,
								   const std::vector<RDGeom::Point3D>& coordinates);
	
		const std::vector<RDGeom::Point3D>& getFeatureCoordinates(FeatureType featureType, const Atom* atom) const
		{
			return featureCoordinates.at(featureType).at(atom);
		}

		Conformer& getConformer()
		{
			return conformer;
		}

		std::map<FeatureType, std::map<const Atom*, std::vector<RDGeom::Point3D>>>& getFeatureCoordinates()
		{
			return featureCoordinates;
		}

		const Conformer& getConformer() const
		{
			return conformer;
		}

		const std::map<FeatureType, std::map<const Atom*, std::vector<RDGeom::Point3D>>>& getFeatureCoordinates() const
		{
			return featureCoordinates;
		}

		std::vector<const Atom*> getFeatureAtoms(FeatureType featureType) const;

        void transformCoordinates(const Eigen::Transform<double, 3, Eigen::Affine> &transform);
	};
}
