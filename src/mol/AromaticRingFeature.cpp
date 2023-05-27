#include "AromaticRingFeature.h"
#include "../gape/SuperpositionMolecule.h"

namespace Gape
{
	thread_local double AromaticRingFeature::normalLength = 3.0;

	std::vector<std::shared_ptr<Feature>> AromaticRingFeature::findAromaticRings(
		const SuperpositionMolecule* superpositionMolecule)
	{
		const auto& mol = superpositionMolecule->getMol();
		const auto ringInfo = mol.getRingInfo();
		int featureSetNumber = 0;
		std::vector<std::shared_ptr<Feature>> features;

		for (const auto ringAtomsIdxs : ringInfo->atomRings())
		{
			const bool aromatic = std::all_of(ringAtomsIdxs.begin(), ringAtomsIdxs.end(), [&mol](const int atomIdx)
			{
				const auto atom = mol.getAtomWithIdx(atomIdx);
				return atom->getIsAromatic();
			});
			const bool planar = std::all_of(ringAtomsIdxs.begin(), ringAtomsIdxs.end(), [&mol](const int atomIdx)
			{
				const auto atom = mol.getAtomWithIdx(atomIdx);
				return atom->getHybridization() == Atom::SP2;
			});

			if (aromatic || planar)
			{
				std::vector<const Atom*> atoms;
				atoms.reserve(ringAtomsIdxs.size());
				std::transform(ringAtomsIdxs.begin(), ringAtomsIdxs.end(), std::back_inserter(atoms),
				               [&mol](const int atomIdx) { return mol.getAtomWithIdx(atomIdx); });
				auto ring = std::make_shared<AromaticRingFeature>(featureSetNumber++, superpositionMolecule, atoms,
				                                                  aromatic, planar);
				features.push_back(std::static_pointer_cast<Feature>(ring));
			}
		}
		return features;
	}

	void AromaticRingFeature::calculateCoordinates(SuperpositionCoordinates& superpositionCoordinates) const
	{
		// find atom to use to store ring feature coordinates.  We need to make sure that it is not in
		// use by another ring.
		const Atom* notInUse = nullptr;
		const Atom* onlyInThisRing = nullptr;
		const auto& mol = molecule->getMol();
		auto atomsInUse = superpositionCoordinates.getFeatureAtoms(FeatureType::AromaticRing);
		for (const auto atom : ringAtoms)
		{
			if (notInUse == nullptr && std::find(atomsInUse.begin(), atomsInUse.end(), atom) == atomsInUse.end())
			{
				notInUse = atom;
			}
			if (mol.getRingInfo()->numAtomRings(atom->getIdx()) == 1)
			{
				onlyInThisRing = atom;
				break;
			}
		}

	}
} // namespace Gape
