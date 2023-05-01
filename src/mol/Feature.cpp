//
// Created by Gareth Jones on 1/19/23.
//

#include <boost/format.hpp>
#include "Feature.h"
#include "MolUtil.h"
#include "../util/Reporter.h"
#include "gape/SuperpositionMolecule.h"

namespace Gape
{
	thread_local double Feature::maximumGaussianScore = 2.0;
	thread_local double Feature::solventVolOk = .0;
	thread_local double Feature::radius = NAN;
	thread_local double Feature::gaussianN = NAN;
	thread_local double Feature::alpha = NAN;

	//const int Feature::N_BUILTIN_FEATURES = 4;
	//const int Feature::MAX_FEATURES = 14;

	void Feature::initialize()
	{
		setRadius(2.0);
	}

	std::string Feature::atomLabel() const
	{
		auto format = boost::format("%s_%d") % featureSetName % (featureSetNumber + 1);
		return format.str();
	}

	std::string Feature::featureLabel(const SuperpositionCoordinates& superpositionCoordinates) const
	{
		auto pharmFeatureGeometry = getPharmFeatureGeometry(superpositionCoordinates);
		auto featureLabel = pharmLabel() + " " + pharmFeatureGeometry->summary();
		/*
		if (pharmacophorePoint)
		{
			featureLabel = featureLabel + " N_MATCHED=" + std::to_string(numberMatched);
		}
		*/
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

	double Feature::solvationPenalty(const RDGeom::Point3D& point,
	                                 const SuperpositionCoordinates& superpositionCoordinates) const
	{
		double penalty = .0;
		const auto& conformer = superpositionCoordinates.getConformer();
		const auto& mol = molecule->getMol();

		for (const auto atom2 : mol.atoms())
		{
			if (atom2->getAtomicNum() == 0)
				continue;
			if (atom2->getIdx() == atom->getIdx())
				continue;
			const double sqrDist = squareDistance(conformer.getAtomPos(atom2->getIdx()), point);
			const double val = -0.5 * solvationAlpha * sqrDist;
			penalty += exp(val);
		}
		return penalty - solventVolOk;
	}

	double Feature::calculateSqrDist(const Feature& otherFeature, const SuperpositionCoordinates& coordinates,
	                                 const SuperpositionCoordinates& otherCoordinates) const
	{
		const auto& point = getFittingPoint(coordinates);
		const auto& otherPoint = otherFeature.getFittingPoint(otherCoordinates);
		const auto squaredDistance = (point - otherPoint).lengthSq();
		return squaredDistance;
	}

	std::string Feature::mappingInfo(const Feature& otherFeature, const SuperpositionCoordinates& coordinates,
	                                 const SuperpositionCoordinates& otherCoordinates) const
	{
		const auto squaredDistance = calculateSqrDist(otherFeature, coordinates, otherCoordinates);
		std::stringstream ss;
		ss << info();
		if (mapped)
		{
			ss << " --> " << otherFeature.info() << " sqrDist " << squaredDistance;
		}
		else
		{
			ss << ":Not Mapped";
		}
		return ss.str();
	}

	double Feature::score(const double sqrDistance)
	{
		return score(sqrDistance, alpha);
	}

	double Feature::score(const double sqrDistance, const double a)
	{
		const double val = -0.5 * a * sqrDistance;
		return exp(val);
	}

	double Feature::getAlpha(const double r)
	{
		return (-2.0 * log(0.5)) / (r * r);
	}
} // Gape
