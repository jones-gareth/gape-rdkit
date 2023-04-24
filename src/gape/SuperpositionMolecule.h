//
// Created by gareth on 10/18/22.
//

#ifndef GAPE_SUPERPOSITIONMOLECULE_H
#define GAPE_SUPERPOSITIONMOLECULE_H

#include <GraphMol/GraphMol.h>
#include <ForceField/MMFF/Params.h>
#include "GapeApp.h"

using namespace RDKit;

namespace ForceFields
{
	class ForceField;
}


namespace RDKit::MMFF
{
	class MMFFMolProperties;
}

namespace Gape
{
	enum class RotatableBondType
	{
		None,
		Flip,
		Full
	};

	class RotatableBond;

	class VdwInfo
	{
	public:
		const unsigned int index0, index1;
		const ForceFields::MMFF::MMFFVdWRijstarEps mmffVdw;

		VdwInfo(unsigned int index0, unsigned int index1, ForceFields::MMFF::MMFFVdWRijstarEps mmffVdw) :
			index0(index0), index1(index1), mmffVdw(mmffVdw)
		{
		}
	};

	class SuperpositionMolecule
	{
	public:
		explicit SuperpositionMolecule(const ROMol& mol, const GapeApp& settings);

		virtual ~SuperpositionMolecule();

		SuperpositionMolecule(const SuperpositionMolecule&) = delete;

		SuperpositionMolecule& operator=(SuperpositionMolecule&) = delete;

		std::string ToMolBlock() const;

		const RWMol& getMol() const { return mol; }

		MMFF::MMFFMolProperties* getMMFFMolProperties() const { return mmffMolProperties; }

		void findFreelyRotatableBonds();

		const std::vector<std::shared_ptr<RotatableBond>>& getRotatableBonds() const { return rotatableBonds; }

		void findPairsToCheck();

		void generate3D();

		void solvate();

		void setConformer(const Conformer& conformer);

		const Conformer& getReferenceConformer() const { return referenceConformer; }

		void findDonorsAndAcceptors();

		const std::map<const Atom*, std::shared_ptr<const HydrogenBondingType>>& getDonors() const { return donors; }

		const std::map<const Atom*, std::shared_ptr<const HydrogenBondingType>>& getAcceptors() const
		{
			return acceptors;
		}

		/**
		 * Check to see if this carbon is in a nitro group.
		 * 
		 * @return
		 */
		bool isNitroOxygen(const Atom& atom) const;

		/**
		 * Checks to see if this oxygen is in a carboxylate group
		 * 
		 * @return
		 */
		bool isCarboxylateOxygen(const Atom& atom) const;

	private:
		RWMol mol;
		MMFF::MMFFMolProperties* mmffMolProperties;
		const GapeApp& settings;
		std::vector<std::shared_ptr<RotatableBond>> rotatableBonds;
		std::vector<VdwInfo> pairsToCheck;
		std::map<const Atom*, std::shared_ptr<const HydrogenBondingType>> donors;
		std::map<const Atom*, std::shared_ptr<const HydrogenBondingType>> acceptors;
		Conformer referenceConformer;

		bool isO2(const Atom& atom) const;

		bool isO3(const Atom& atom) const;

		bool isAmideBond(const Bond& bond) const;

		bool isNpl3Atom(const Atom& atom) const;

		static bool isTerminalBond(const Bond& bond);

		bool isArginineCarbon(const Atom& atom) const;

		static bool isSp2Carbon(const Atom& atom);

		bool atomIsInRing(const Atom& atom) const;

		bool isCOOHCarbon(const Atom& atom, Atom*& o2Atom, Atom*& o3Atom) const;

		RotatableBondType isRotatableBond(const Bond& bond, bool& canFlatten) const;

		/**
		 * Check to see if this carbon is in a nitro group.
		 * 
		 * @return
		 */
		bool isNitroNitrogen(const Atom& atom) const;

		/**
		 * Checks to see if this carbon is a carboxylate
		 * 
		 * @return
		 */
		bool isCarboxylateCarbon(const Atom& atom) const;
	};
} // namespace GAPE

#endif //GAPE_SUPERPOSITIONMOLECULE_H
