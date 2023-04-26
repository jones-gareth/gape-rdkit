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

		/**
		 * \brief Number lone pairs
		 */
		int numberLonePairs;

		/**
		 * \brief Hydrogen bonding type of acceptor
		 */
		std::shared_ptr<const HydrogenBondingType> hydrogenBondingType;

		/**
		 * \brief Reference coordinates of lone pairs
		 */
		std::vector<RDGeom::Point3D> lonePairs;

		static const double lonePairLength;

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
		static void addThreePairsToTetrahedral(const RDGeom::Point3D& origin, const RDGeom::Point3D& atom1,
		                                       RDGeom::Point3D& lonePair1, RDGeom::Point3D& lonePair2,
		                                       RDGeom::Point3D& lonePair3,
		                                       double lonePairLength);
		/**
		 * Acceptor is at origin and is bonded to atom1 and atom2. Two lone pairs
		 * are added to complete the sp3 geometry.
		 * 
		 * @param origin
		 * @param atom1
		 * @param atom2
		 * @param lp1
		 * @param lp2
		 * @param lonePairLength
		 */
		static void addTwoPairsToTetrahedral(const RDGeom::Point3D& origin, const RDGeom::Point3D& atom1,
		                                     const RDGeom::Point3D& atom2, RDGeom::Point3D& lonePair1,
		                                     RDGeom::Point3D& lonePair2, double lonePairLength);


		/**
		 * Add one lone pair to an acceptor bonded to three atoms in sp3 geometry.
		 * The acceptor is at co-ordinates atom and three atoms (atom2, atom3 and
		 * atom4) are bonded to it. If C is the centroid of atoms 2, 3 and 4 and C'
		 * is reflection through the origin, O, then the lone pair is added on OC'
		 * 
		 * @param atom
		 * @param atom2
		 * @param atom3
		 * @param atom4
		 * @param lp
		 * @param lonePairLength
		 */
		static void addOnePairToTetrahedral(const RDGeom::Point3D& origin, const RDGeom::Point3D& atom1,
		                                    const RDGeom::Point3D& atom2, const RDGeom::Point3D& atom3,
		                                    RDGeom::Point3D& lonePair1, double lonePairLength);

		/**
		 * Add one lone pair to the sp2 acceptor with co-ordinates origin. The
		 * acceptor is bonded to atoms atom1 and atom2.
		 * 
		 * @param origin
		 * @param atom1
		 * @param atom2
		 * @param lp
		 * @param lonePairLength
		 * @return
		 */
		static bool addOnePairToTrigonal(const RDGeom::Point3D& origin, const RDGeom::Point3D& atom1,
		                                 const RDGeom::Point3D& atom2, RDGeom::Point3D& lonePair1,
		                                 double lonePairLength);

		/**
		 * Adds two pair to the acceptor with sp2 geometry at co-ordinates <atom1>.
		 * The acceptor is bonded to <origin> and <origin> is connected to atom2.
		 * 
		 * @param atom1
		 * @param origin
		 * @param atom2
		 * @param lp1
		 * @param lp2
		 * @param lonePairLength
		 * @return
		 */
		static bool addTwoPairsToTrigonal(const RDGeom::Point3D& atom1, const RDGeom::Point3D& origin,
		                                  const RDGeom::Point3D& atom2, RDGeom::Point3D& lonePair1,
		                                  RDGeom::Point3D& lonePair2, double lonePairLength);

		/**
		 * Adds two lone pairs to the sp2 acceptor atom1 which is bonded to the atom
		 * at co-ordinates origin. The plane in which the lone pairs are added is
		 * undefined
		 * 
		 * alpha argument is double the angle between the lone pairs, so we can use
		 * this transform for the cone feature
		 * 
		 * @param atom1
		 * @param origin
		 * @param lp1
		 * @param lp2
		 * @param alpha
		 * @param lonePairLength
		 */
		static void addTwoPairsRandomlyToTrigonal(const RDGeom::Point3D& atom1, const RDGeom::Point3D& origin,
		                                          RDGeom::Point3D& lonePair1, RDGeom::Point3D& lonePair2, double alpha,
		                                          double lonePairLength);

		/**
		 * Returns the true number of lone pairs that an acceptor atom has. Differs
		 * from countLonePairs in that an acceptor with no geometry or cone geometry
		 * will return the true number of lone pairs (probably 3).
		 * 
		 * @param a
		 * @return
		 */
		int countTrueLonePairs() const;

		/**
		 * Returns a count of lone pairs that need to be added to an atom. Take
		 * account of a reduced representation whereby an acceptor that accepts in a
		 * cone is represented by one (rather that 3) lone pairs
		 * 
		 * @param a
		 * @return
		 */
		int countLonePairs() const;

		/**
		 * Adds one lone pair to an sp3 acceptor.
		 * 
		 */
		void add1LpToSp3();

		/**
		 * Adds one lone pair to an sp2 atom.
		 * 
		 */
		void add1LpToSp2();

		/**
		 * Adds one lone pair to an sp1 atom.
		 * 
		 */
		void addLpToSp1();

		/**
		 * Adds two lone pairs to an sp3 atom.
		 * 
		 */
		void add2LpToSp3();

		/**
		 * Adds two lone pairs to an sp2 atom.
		 * 
		 */
		void add2LpToSp2();

		/**
		 * Adds two lone pairs to the carboxylic oxygen O.co2 atom. Also used to add
		 * lone pairs to nitro oxygen.
		 * 
		 */
		void addLpToOco2();


		/**
		 * Adds three lone pairs to an sp3 atom.
		 * 
		 */
		void add3LpToSp3();

		/**
		 * \brief Add lone pairs to acceptor based on atom hybridization
		 */
		void addAtomLonePairs();

		/**
		 * \brief Add one lone pair to atom regardless of atom hybridization
		 */
		void addOnePairToAtom();

		/**
		 * \brief Add lone pairs to atom based on acceptor geometry
		 */
		void addLonePairs();

		/**
		 * \brief 
		 * \return validate lone pair lengths are correct
		 */
		bool checkLonePairLengths();

	public:
		AcceptorAtom(const SuperpositionMolecule* spMol, const Atom* atom);

		const Atom* getAtom() const { return atom; }
	};
}
