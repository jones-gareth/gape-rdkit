//
// Created by gareth on 10/18/22.
//

#pragma once

#include <memory>
#include <vector>

#include "SuperpositionMolecule.h"

namespace Gape
{
	class SuperpositionGa
	{
		const static int MAX_MOLS;
		SuperpositionMolecule *baseMolecule, *fittingMolecule;

		std::vector<SuperpositionMolPtr> molecules;
		const GapeSettings& settings;

		SuperpositionGa(const std::vector<std::shared_ptr<SuperpositionMolecule>>& m,
		                const GapeSettings& s);

		void setupMolecules();
		void findBaseMolecule();

	private:
		SuperpositionMolecule* leastFlexibleMolecule() const;
		SuperpositionMolecule* mostFlexibleMolecule() const;
		SuperpositionMolecule* mostFeaturedMolecule() const;
		SuperpositionMolecule* mostFeaturedRigidMolecule() const;
		SuperpositionMolecule* mostActiveRigidMolecule() const;
		SuperpositionMolecule* mostActiveMolecule() const;

		size_t numberMoleculesMapped;
		std::vector<int> integerEntryPoints;
		int integerStringLength;
		int baseMoleculeNumber = -1;
		int fittingMoleculeNumber = -1;

	};
} // GapeApp
