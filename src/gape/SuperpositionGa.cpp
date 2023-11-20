//
// Created by gareth on 10/18/22.
//

#include "SuperpositionGa.h"

#include "util/Reporter.h"

namespace Gape
{
	const int SuperpositionGa::MAX_MOLS = 100;

	SuperpositionGa::SuperpositionGa(const std::vector<std::shared_ptr<SuperpositionMolecule>>& m,
	                                 const GapeSettings& s): molecules(m), settings(s)
	{
		setupMolecules();
	}

	void SuperpositionGa::setupMolecules()
	{
		if (molecules.size() > MAX_MOLS)
		{
			molecules.resize(MAX_MOLS);
		}
	}

	void SuperpositionGa::findBaseMolecule()
	{
		baseMolecule = molecules[0].get();
		fittingMolecule = molecules[0].get();

		bool rigidBase = false, multipleRigid = false;
		for (const auto& molecule : molecules)
		{
			if (molecule->isRigid())
			{
				if (rigidBase)
				{
					multipleRigid = true;
					baseMolecule = molecules[0].get();
					fittingMolecule = molecules[0].get();
					rigidBase = false;
					break;
				}
				REPORT(Reporter::INFO) << "Base molecule is rigid";
				baseMolecule = fittingMolecule = molecule.get();
				rigidBase = true;
			}
		}
		// I thought for feature clustering it would be best to use the least
		// flexible molecule as the fitting molecule, but this doesn't seem to
		// work in practice.
		constexpr bool featureClusteringUseLeastFlexible = false;
		const auto baseMoleculeSelection = settings.getGapeParameters().baseMoleculeSelection;

		// Base molecule selection
		if (rigidBase)
		{
		}
		else if (multipleRigid)
		{
			// in the case of multiple rigid use max features when min rotatable
			// bonds is selected
			if (baseMoleculeSelection == BaseMoleculeSelection::minRotatableBonds ||
				baseMoleculeSelection == BaseMoleculeSelection::maxFeatures)
			{
				REPORT(Reporter::INFO) << "Base molecule selected on maximum features";
				baseMolecule = mostFeaturedRigidMolecule();
			}
			else if (baseMoleculeSelection == BaseMoleculeSelection::maxActivity)
			{
				REPORT(Reporter::INFO) << "Base molecule selected on maximum activity";
				if (!settings.getGapeParameters().useActivities)
				{
					const std::string message(
						"Base molecule selection method is maximum activity, but use_activities is false");
					REPORT(Reporter::FATAL) << message;
					throw runtime_error(message);
				}
				baseMolecule = mostActiveRigidMolecule();
			}
			else
			{
				const std::string message("Unknown selection strategy ");
				REPORT(Reporter::FATAL) << message;
				throw runtime_error(message);
			}
		}

		else if (featureClusteringUseLeastFlexible)
		{
			REPORT(Reporter::INFO) << "Using feature clustering - base molecule is least flexible";
			baseMolecule = leastFlexibleMolecule();
		}

		else
		{
			if (baseMoleculeSelection == BaseMoleculeSelection::minRotatableBonds)
			{
				REPORT(Reporter::INFO) << "Base molecule selected on minimum rotatable bonds";
				baseMolecule = leastFlexibleMolecule();
			}
			else if (baseMoleculeSelection == BaseMoleculeSelection::maxFeatures)
			{
				REPORT(Reporter::INFO) << "Base molecule selected on maximum features";
				baseMolecule = mostFeaturedMolecule();
			}
			else if (baseMoleculeSelection == BaseMoleculeSelection::maxActivity)
			{
				REPORT(Reporter::INFO) << "Base molecule selected on maximum activity";
				if (!settings.getGapeParameters().useActivities)
				{
					const std::string message(
						"Base molecule selection method is maximum activity, but use_activities is false");
					REPORT(Reporter::FATAL) << message;
					throw runtime_error(message);
				}
				baseMolecule = mostActiveMolecule();
			}
			else
				throw std::runtime_error("Unknown selection strategy");
		}

		const auto fittingMoleculeSelection = settings.getGapeParameters().fittingMoleculeSelection;
		if (rigidBase)
		{
			;
		}

		else if (multipleRigid)
		{
			if (fittingMoleculeSelection == FittingMoleculeSelection::maxFeatures ||
				fittingMoleculeSelection == FittingMoleculeSelection::minRotatableBonds)
			{
				REPORT(Reporter::INFO) << "Fitting molecule selected on maximum features";
				fittingMolecule = mostFeaturedRigidMolecule();
			}
			else if (fittingMoleculeSelection == FittingMoleculeSelection::baseMolecule)
			{
				REPORT(Reporter::INFO) << "Fitting molecule is base molecule";
				fittingMolecule = baseMolecule;
			}
			else
				throw runtime_error("Unknown selection strategy");
		}

		else if (featureClusteringUseLeastFlexible)
		{
			REPORT(Reporter::INFO) << "Using feature clustering - fitting molecule is least flexible";
			fittingMolecule = leastFlexibleMolecule();
		}

		else
		{
			if (fittingMoleculeSelection == FittingMoleculeSelection::minRotatableBonds)
			{
				REPORT(Reporter::INFO) << "Fitting molecule selected on minimum rotatable bonds";
				fittingMolecule = leastFlexibleMolecule();
			}
			else if (fittingMoleculeSelection == FittingMoleculeSelection::maxFeatures)
			{
				REPORT(Reporter::INFO) << "Fitting molecule selected on maximum features";
				fittingMolecule = mostFeaturedMolecule();
			}
			else if (fittingMoleculeSelection == FittingMoleculeSelection::baseMolecule)
			{
				REPORT(Reporter::INFO) << "Fitting molecule is base molecule";
				fittingMolecule = baseMolecule;
			}
			else
				throw std::runtime_error("Unknown selection strategy");
		}
		REPORT(Reporter::INFO) << "Base molecule is " << baseMolecule->getName();


		const std::function isMapped = [this](const SuperpositionMolPtr& mol)
		{
			return !(mol->isFixed() || mol.get() == fittingMolecule);
		};
		const auto mappedMolecules = filterListToNewList(molecules, isMapped);
		numberMoleculesMapped = mappedMolecules.size();
		const std::vector<int> iep(numberMoleculesMapped);
		integerEntryPoints = iep;
		integerStringLength = static_cast<int>(fittingMolecule->numberMappingFeatures() * numberMoleculesMapped);

		int no = 0;
		for (size_t i = 0; i < molecules.size(); ++i)
		{
			const auto& molecule = molecules[i];
			if (molecule.get() == baseMolecule)
			{
				baseMoleculeNumber = static_cast<int>(i);
			}
			if (molecule.get() == fittingMolecule)
			{
				fittingMoleculeNumber = static_cast<int>(i);
				continue;
			}
			if (molecule->isFixed())
				continue;
			integerEntryPoints[no] = no * static_cast<int>(molecule->numberMappingFeatures());
			++no;
		}
		/*
				integerStringRanges = new int[integerStringLength];
				int pos = 0;
				Map<FeatureType, FeatureMapping> fittingFeatures = fittingMolecule
						.getFeatureMappings();
				for (GaMolecule molecule : molecules) {
					if (molecule == fittingMolecule)
						continue;
					if (molecule.isFixed())
						continue;
					for (Entry<FeatureType, FeatureMapping> entry : fittingFeatures.entrySet()) {
						if (!entry.getValue().isMapping())
							continue;
						FeatureMapping otherFeatures = molecule
								.getFeatureMappings(entry.getKey());
						for (int k = 0; k < entry.getValue().getNFeatures(); k++) {
							integerStringRanges[pos] = otherFeatures.getNFeatures();
							pos++;
						}
					}
				}
		
				setWeights();
			*/
	}


	SuperpositionMolecule* SuperpositionGa::leastFlexibleMolecule() const
	{
		const std::function score = [](const SuperpositionMolPtr& mol) { return -mol->getRotatableBonds().size(); };
		const auto& leastFlexible = findMaxBy(molecules, score);
		return leastFlexible.get();
	}

	SuperpositionMolecule* SuperpositionGa::mostFeaturedMolecule() const
	{
		const std::function score = [](const SuperpositionMolPtr& mol) { return mol->numberFeatures(); };
		const auto& mostFeatured = findMaxBy(molecules, score);
		return mostFeatured.get();
	}

	SuperpositionMolecule* SuperpositionGa::mostFlexibleMolecule() const
	{
		const std::function score = [](const SuperpositionMolPtr& mol) { return mol->getRotatableBonds().size(); };
		const auto& mostFlexible = findMaxBy(molecules, score);
		return mostFlexible.get();
	}

	SuperpositionMolecule* SuperpositionGa::mostFeaturedRigidMolecule() const
	{
		const std::function filter = [](const SuperpositionMolPtr& mol) { return mol->isRigid(); };
		const auto rigidMolecules = filterListToNewList(molecules, filter);
		const std::function score = [](const SuperpositionMolPtr& mol) { return mol->numberFeatures(); };
		const auto& mostFeatured = findMaxBy(rigidMolecules, score);
		return mostFeatured.get();
	}

	Gape::SuperpositionMolecule* Gape::SuperpositionGa::mostActiveRigidMolecule() const
	{
		const std::function filter = [](const SuperpositionMolPtr& mol) { return mol->isRigid(); };
		const auto rigidMolecules = filterListToNewList(molecules, filter);
		const std::function score = [](const SuperpositionMolPtr& mol) { return mol->activity(); };
		const auto& mostActive = findMaxBy(rigidMolecules, score);
		return mostActive.get();
	}

	Gape::SuperpositionMolecule* Gape::SuperpositionGa::mostActiveMolecule() const
	{
		const std::function score = [](const SuperpositionMolPtr& mol) { return mol->activity(); };
		const auto& mostActive = findMaxBy(molecules, score);
		return mostActive.get();
	}
} // GapeApp
