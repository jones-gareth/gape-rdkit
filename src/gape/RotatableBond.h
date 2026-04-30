//
// Created by Gareth Jones on 11/25/2022.
//

#ifndef GAPE_ROTATABLEBOND_H
#define GAPE_ROTATABLEBOND_H

#include <GraphMol/GraphMol.h>
#include "SuperpositionMolecule.h"
#include <ForceField/MMFF/Contribs.h>
#include <Geometry/Transform3D.h>
#include "TorsionInfo.h"

namespace Gape
{
	using namespace RDKit;

	class SuperpositionMolecule;
	class CornerRotation;

	class RotatableBond
	{
	private:
		const Atom *atom1, *atom2;
		RotatableBondType rotatableBondType;
		const SuperpositionMolecule* molecule;
		std::vector<const Atom*> atom1List, atom2List;
		std::vector<TorsionInfo> torsions;

		void setTorsionAngles();

		friend SuperpositionMolecule;
		friend CornerRotation;

	public:
		RotatableBond(RotatableBondType rotatableBondType, const Bond* bond,
		              const SuperpositionMolecule* superpositionMolecule);

		bool isSeparatedByBond(const Atom* a1, const Atom* a2) const;

		RotatableBondType getRotatableBondType() const { return rotatableBondType; }

		void rotateBond(double angle, Conformer& conf) const;

		void rotateBond(double angle, SuperpositionCoordinates& superpositionCoordinates, RDGeom::Transform3D &rot) const;

		void rotateBond(double angle, SuperpositionCoordinates& superpositionCoordinates) const;

        double rotatableBondEnergy(const Conformer &conformer) const;

		const std::vector<TorsionInfo>& getTorsions() const { return torsions; }

		const std::vector<const Atom*>& getAtom1List() const { return atom1List; }

		const std::vector<const Atom*>& getAtom2List() const { return atom2List; }

		~RotatableBond();

		void flattenBond(Conformer &conformer) const;
	};
}

#endif //GAPE_ROTATABLEBOND_H
