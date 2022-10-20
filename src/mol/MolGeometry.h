/*
 * MolGeometry.h
 *
 * Routines for determining co-ordinates for adding hydrogens or lone pairs to atoms
 *
 *  Created on: Dec 16, 2015
 *      Author: Gareth Jones
 */

#ifndef SRC_MOL_MOLGEOMETRY_H_
#define SRC_MOL_MOLGEOMETRY_H_

#include "Atom.h"
#include "Molecule.h"

/**
 * Functions to perform molecular geometry operations such as determining cooordinates
 * for lone pairs or hydrogens.
 *
 * !!!! Most of these routines have not been tested yet.
 */
namespace GarethMol {

const double DEFAULT_DISTANCE = 1.0;

/**
 * Returns the coordinates of a single point which is an linear extension of the
 * vector from atom1 to atom2 at the specified distance from atom1.
 *
 * @param molecule
 * @param atom1
 * @param atom2
 * @param distance
 * @return
 */
CoordVector generateOnePointForLinear(const Molecule & molecule,
		const Atom & atom1, const Atom & atom2, const double distance =
				DEFAULT_DISTANCE);

/**
 * Returns the three remaining tetrahedral points (at the specified distance from
 * atom1) for the tetrahedral system centered on atom1.  Atom2 defines the forth
 * vertex in the system.
 *
 * @param molecule
 * @param atom1
 * @param atom2
 * @param distance
 * @return
 */
std::vector<CoordVector> generateThreePointsForTetrahedral(const Molecule & molecule,
		const Atom & atom1, const Atom & atom2, const double distance =
				DEFAULT_DISTANCE);
/**
 * Returns the two remaining tetrahedral points (at the specified distance from
 * atom1) for the tetrahedral system centered on atom1.  Atom2 and atom3 define
 * the third and forth vertices in the system.
 *
 * @param molecule
 * @param atom1
 * @param atom2
 * @param distance
 * @return
 */
std::vector<CoordVector> generateTwoPointsForTetrahedral(const Molecule & molecule,
		const Atom & atom1, const Atom & atom2, const Atom & atom3,
		const double distance = DEFAULT_DISTANCE);

/**
 * Returns the remaining tetrahedral point (at the specified distance from
 * atom1) for the tetrahedral system centered on atom1.  Atom2, atom3 and atom4 define
 * the other vertices in the system.
 *
 * @param molecule
 * @param atom1
 * @param atom2
 * @param atom3
 * @param distance
 * @return
 */
CoordVector generateOnePointForTetrahedral(const Molecule & molecule,
		const Atom & atom1, const Atom & atom2, const Atom & atom3,
		const Atom & atom4, const double distance = DEFAULT_DISTANCE);

/**
 * Returns the remaining trigonal point (at the specified distance from
 * atom1) for the trigonal system centered on atom1.  Atom2 and atom3 define
 * the other vertices in the system.
 *
 * @param molecule
 * @param atom1
 * @param atom2
 * @param atom3
 * @param distance
 * @return
 */
CoordVector generateOnePointForTrigonal(const Molecule & molecule,
		const Atom & atom1, const Atom & atom2, const Atom & atom3,
		const double distance = DEFAULT_DISTANCE);

/**
 * Returns the remaining two trigonal points for atom1.  Atom 1 is bonded to atom2.
 * Atom2 is further bonded to atom 3 which defines the plane of the trigonal system.
 *
 * @param molecule
 * @param atom1
 * @param atom2
 * @param atom3
 * @param distance
 * @return
 */
std::vector<CoordVector> generateTwoPointsForTrigonal(const Molecule & molecule,
		const Atom & atom1, const Atom & atom2, const Atom & atom3,
		const double distance = DEFAULT_DISTANCE);

/**
 * Returns the remaining two trigonal points for atom1.  Atom 1 is bonded to atom2,
 * but the plane of the trigonal system is undefined.
 *
 * @param molecule
 * @param atom1
 * @param atom2
 * @param distance
 * @return
 */
std::vector<CoordVector> generateTwoPointsForIncompleteTrigonal(
		const Molecule & molecule, const Atom & atom1, const Atom & atom2,
		const double distance = DEFAULT_DISTANCE);

double outOfPlaneHeight(const Molecule & molecule, const Atom & atom,
		const Atom & neighbour1, const Atom & neighbour2, const Atom &neighbour3);

}
/* namespace GarethMol */

#endif /* SRC_MOL_MOLGEOMETRY_H_ */
