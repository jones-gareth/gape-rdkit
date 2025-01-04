#include "AromaticRingFeature.h"
#include "gape/SuperpositionMolecule.h"
#include "util/Reporter.h"


namespace Gape
{
	namespace Detail
	{
		std::vector<int> atomsToAtomIds(const std::vector<const Atom*>& atoms)
		{
			std::vector<int> atomIds;
			atomIds.reserve(atoms.size());
			std::transform(atoms.begin(), atoms.end(), std::back_inserter(atomIds),
			               [](const Atom* atom) { return atom->getIdx(); });
			std::sort(atomIds.begin(), atomIds.end());
			return atomIds;
		}
	}

	thread_local double AromaticRingFeature::normalLength = 3.0;

	std::vector<std::shared_ptr<Feature>> AromaticRingFeature::findAromaticRings(
		const SuperpositionMolecule* superpositionMolecule)
	{
		const auto& mol = superpositionMolecule->getMol();
		const auto ringInfo = mol.getRingInfo();
		int featureSetNumber = 0;
		std::vector<std::shared_ptr<Feature>> features;

		std::vector<const Atom*> atomsInUse;
		for (const auto& ringAtomsIdxs : ringInfo->atomRings())
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
				                                                  aromatic, planar, atomsInUse);
				atomsInUse.push_back(ring->atom);
				features.push_back(std::static_pointer_cast<Feature>(ring));
			}
		}
		return features;
	}

	AromaticRingFeature::AromaticRingFeature(const int featureSetNumber, const SuperpositionMolecule* spMol,
	                                         const std::vector<const Atom*>& ringAtoms, const bool aromatic,
	                                         const bool planar, const std::vector<const Atom*>& atomsInUse):
		AromaticRingFeature(featureSetNumber)
	{
		this->ringAtoms = ringAtoms;
		this->aromatic = aromatic;
		this->planar = planar;
		molecule = spMol;

		// find atom to use to store ring feature coordinates.  We need to make sure that it is not in
		// use by another ring.
		const Atom* notInUse = nullptr;
		const Atom* onlyInThisRing = nullptr;
		const auto& mol = molecule->getMol();
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
		atom = onlyInThisRing != nullptr ? onlyInThisRing : notInUse;
		if (atom == nullptr)
		{
			throw std::domain_error("Unable to find representative atom for planar ring");
		}
	}

	RDGeom::Point3D AromaticRingFeature::ringCenter(const SuperpositionCoordinates& superpositionCoordinates) const
	{
		const auto& conformer = superpositionCoordinates.getConformer();
		RDGeom::Point3D center;
		for (const auto atom : ringAtoms)
		{
			center += conformer.getAtomPos(atom->getIdx());
		}
		center /= static_cast<int>(ringAtoms.size());
		return center;
	}

	std::vector<RDGeom::Point3D> AromaticRingFeature::ringNormals(
		const SuperpositionCoordinates& superpositionCoordinates) const
	{
		const auto& conformer = superpositionCoordinates.getConformer();
		const auto center = ringCenter(superpositionCoordinates);
		RDGeom::Point3D normal;
		for (size_t i = 0; i < ringAtoms.size(); i++)
		{
			const auto& point1 = conformer.getAtomPos(ringAtoms[i]->getIdx());
			const auto atom2 = i == ringAtoms.size() - 1 ? ringAtoms[0] : ringAtoms[i + 1];
			const auto& point2 = conformer.getAtomPos(atom2->getIdx());
			const auto a1 = center - point1;
			const auto a2 = center - point2;
			const auto vp = a1.crossProduct(a2);
			normal += vp;
		}

		normal.normalize();
		normal *= normalLength;

		const auto normal1 = center + normal;
		const auto normal2 = center - normal;

		std::vector<RDGeom::Point3D> centerAndNormals(2);
		centerAndNormals.push_back(center);
		centerAndNormals.push_back(normal1);
		centerAndNormals.push_back(normal2);
		return centerAndNormals;
	}

	void AromaticRingFeature::calculateCoordinates(SuperpositionCoordinates& superpositionCoordinates) const
	{
		const auto centerAndNormals = ringNormals(superpositionCoordinates);
		superpositionCoordinates.addFeatureCoordinates(FeatureType::AromaticRing, atom, centerAndNormals);
	}

	FeatureScore AromaticRingFeature::score(const Feature& otherFeature, const SuperpositionCoordinates& coordinates,
	                                  const SuperpositionCoordinates& otherCoordinates) const
	{
		const auto& other = dynamic_cast<const AromaticRingFeature&>(otherFeature);
		const auto& centerAndNormals = coordinates.getFeatureCoordinates(FeatureType::AromaticRing, atom);
		const auto& otherCenterAndNormals = coordinates.getFeatureCoordinates(FeatureType::AromaticRing, other.atom);

		// Center score
		const auto& center = centerAndNormals[0];
		const auto& otherCenter = otherCenterAndNormals[1];
		const auto centerSqrDistance = (center - otherCenter).lengthSq();
		const auto vol = Feature::score(centerSqrDistance);

		// Normal score
		double nVol = .0;
		for (int i = 1; i < 3; i++)
		{
			const auto& normal = centerAndNormals[i];
			double v = .0;
			for (int j = 1; j < 3; j++)
			{
				const auto& otherNormal = otherCenterAndNormals[j];
				const auto sqrDistance = (normal - otherNormal).lengthSq();
				const auto t = Feature::score(sqrDistance);
				if (t > v)
				{
					v = t;
				}
			}
			nVol += v;
		}
		nVol = nVol / 2;

		// the score is the lesser of center and normal scores
		double score = (vol < nVol) ? vol : nVol;
		REPORT(Reporter::DEBUG) << " Score: " << score << " vol " << vol << " nVol " << nVol;

		if (score > maximumGaussianScore)
		{
			score = maximumGaussianScore;
		}

		REPORT(Reporter::DEBUG) << info() << " " << other.info() << " score " << score;
		const FeatureScore result(score, score);
		return result;
	}

	std::string AromaticRingFeature::pharmLabel() const
	{
		const auto atomIds = Detail::atomsToAtomIds(ringAtoms);
		stringstream ss;
		ss << featureSetName << "ATOMS";
		for (const auto idx : atomIds)
		{
			ss << "_" << idx;
		}
		return ss.str();
	}

	std::string AromaticRingFeature::info() const
	{
		const auto atomIds = Detail::atomsToAtomIds(ringAtoms);
		std::stringstream ss;
		auto first = true;
		for (const int idx : atomIds)
		{
			if (!first)
			{
				ss << ",";
			}
			first = false;
			ss << idx;
		}
		ss << "]";
		return ss.str();
	}

	std::unique_ptr<PharmFeatureGeometry> AromaticRingFeature::getPharmFeatureGeometry(
		const SuperpositionCoordinates& superpositionCoordinates) const
	{
		const auto& centerAndNormals = superpositionCoordinates.getFeatureCoordinates(FeatureType::AromaticRing, atom);
		return std::make_unique<VectorPharmFeatureGeometry>(centerAndNormals[1], centerAndNormals[2]);
	}

	const RDGeom::Point3D& AromaticRingFeature::getFittingPoint(
		const SuperpositionCoordinates& superpositionCoordinates) const
	{
		const auto& centerAndNormals = superpositionCoordinates.getFeatureCoordinates(FeatureType::AromaticRing, atom);
		return centerAndNormals[0];
	}
} // namespace Gape
