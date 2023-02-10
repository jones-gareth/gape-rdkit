//
// Created by Gareth Jones on 1/19/23.
//

#include <boost/format.hpp>
#include "Feature.h"
#include "MolUtil.h"
#include "../util/Reporter.h"

namespace Gape
{
	thread_local double Feature::maximumGaussianScore = 2.0;
	thread_local double Feature::solventVolOk = .0;

	//const int Feature::N_BUILTIN_FEATURES = 4;
	//const int Feature::MAX_FEATURES = 14;

	void Feature::initialize()
	{
		setRadius(2.0);
	}

	std::string Feature::atomLabel() const
	{
		int no = featureSetNumber + 1;
		auto format = boost::format("%s_%d") % featureSetName % featureSetNumber;
		return format.str();
	}

	std::string Feature::featureLabel() const
	{
		auto featureLabel = pharmLabel() + " " + pharmFeatureGeometry->summary();
		if (pharmacophorePoint)
		{
			featureLabel = featureLabel + " N_MATCHED=" + std::to_string(numberMatched);
		}
		return featureLabel;
	}

	void Feature::setRadius(const double r)
	{
		alpha = getAlpha(r);
		gaussianN = pow((2.0 * alpha) / M_PI, 0.75);
		radius = r;
		REPORT(Reporter::DETAIL) << "Setting alpha " << alpha << " gaussianN " << gaussianN
                                 << " radius " << r;
	}

	static double getAlpha(const double r)
	{
		return (-2.0 * log(0.5)) / (r * r);
	}

	double Feature::solvationAlpha = Feature::getAlpha(1.0);

	void Feature::setSolvationAlpha(const double r)
	{
		solvationAlpha = getAlpha(r);
	}

	double Feature::solvationPenalty(const RDGeom::Point3D& point, const ROMol& mol, const Conformer& conformer,
	                                 const Atom& atom) const
	{
		double penalty = .0;

		for (const auto atom2 : mol.atoms())
		{
			if (atom2->getAtomicNum() == 0)
				continue;
			if (atom2->getIdx() == atom.getIdx())
				continue;
			const double sqrDist = squareDistance(conformer.getAtomPos(atom2->getIdx()), point);
			const double val = -0.5 * solvationAlpha * sqrDist;
			penalty += exp(val);
		}
		return penalty - solventVolOk;
	}

	const RDGeom::Point3D& Feature::calculateCoordinate(const Conformer& conformer)
	{
		coordinate = conformer.getAtomPos(atom->getIdx());
		return coordinate;
	}

	const RDGeom::Point3D& Feature::getSavedCoordinate() const
	{
		return coordinate;
	}

	double Feature::calculateSqrDist()
	{
		squaredDistance = (coordinate - mappedFeature->coordinate).lengthSq();
		return squaredDistance;
	}

	std::string Feature::mappingInfo() const
	{
		std::stringstream ss;
		ss << info();
		if (mapped)
		{
			ss << " --> " << mappedFeature->info() << " sqrDist " << squaredDistance;
		}
		else
		{
			ss << ":Not Mapped";
		}
		return ss.str();
	}
} // Gape
