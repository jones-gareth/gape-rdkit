#pragma once


#include <GraphMol/GraphMol.h>

#include "../gape/SuperpositionMolecule.h"

using namespace RDKit;

namespace Gape
{
	class AcceptorAtom
	{
		/**
		 * Molecule containing this feature
		 */
		const SuperpositionMolecule* molecule;

		/**
		 * Acceptor atom
		 */
		const Atom* atom;

		int numberLonePairs;

		/**
		 * Adds one LP to atom1 on a line that runs from atom2 through atom1.
		 */
		static void addOnePairToLinear(const RDGeom::Point3D& atom1, const RDGeom::Point3D& atom2,
		                               RDGeom::Point3D& lonePair,
		                               double lonePairLength);

		/**
		 * origin is the acceptor atom and atom1 is bonded to it. Three lone pairs are
		 * added to complete sp3 geometry.
		 */
		static void addThreePairsToTetrahedral(const RDGeom::Point3D& origin, const RDGeom::Point3D& atom2,
		                                       RDGeom::Point3D& lonePair1, RDGeom::Point3D& lonePair2,
		                                       RDGeom::Point3D& lonePair3,
		                                       double lonePairLength);

	public:
		AcceptorAtom(const SuperpositionMolecule* spMol, const Atom* atom);
	};
}
