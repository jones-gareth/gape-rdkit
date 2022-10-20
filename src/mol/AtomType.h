/*
 * AtomType.h
 *
 *  Created on: Jun 21, 2013
 *      Author: Gareth Jones
 *
 *  Defined available atom types
 */

#include <string>
#include <cmath>
#include <vector>
#include <map>
#include <memory>

#ifndef ATOMTYPE_H_
#define ATOMTYPE_H_

namespace GarethMol {

class AtomType {

public:

	/**
	 * Enumerated list of atom type identifiers
	 */
	enum class AtomTypeId {
		ATM_NONE,

		C3, C2, C1, CAR, CCAT,

		N3, N2, N1, NAR, NAM, NPL3, N4,

		O3, O2, OCO2, OAR,

		S3, S2, SO, SO2, P3, H, CL, BR, B,

		I, SI, NA, K, F, CA, LI, AL, ANY, MG, ZN, FE, MN,

		ATOM_USER1, ATOM_USER2, ATOM_USER3, ATOM_USER4, ATOM_USER5,

		ATOM_USER6, ATOM_USER7, ATOM_USER8, ATOM_USER9, ATOM_USER10,

		HET, HAL, HEV,

		C, N, P, S, OANY,

		X, Y, Z, DU, WILD, DUC, LP,

		CU, SE, PB, HG, SB, BI, TI, AS, NI, BA,

		AG, SN, CR, CD, CO, AU, CE, BE, TH, V, ZR,

		PT, TE, RA, GE, GA, TL, CS, IN, SM, PD

	// Most of the metals are there to make sure we can read first 20K smiles in the NCI database
	};

	/**
	 * Geometry types
	 */
	enum class Geometry {
		LIN, TRI, TET, CO, NONE
	};

	virtual ~AtomType();

	/**sdChiral = true;
	 * Return the atom type with a particular name
	 *
	 * @param name
	 * @return
	 */
	static const AtomType & typeFromName(const std::string & name);

	/**
	 * Return the string name for a particular atom type id
	 *
	 * @param type
	 * @return
	 */
	static const std::string & nameFromTypeId(const AtomTypeId & type);

	/**
	 * Return the Atom type for a particular atom type id
	 *
	 * @param type
	 * @return
	 */
	static const AtomType & typeFromTypeId(const AtomTypeId & type);

	/**
	 * Return the atom type for a particular atomic number
	 * @param atomicNumber
	 * @return
	 */
	static const AtomType & typeFromAtomicNumber(const int atomicNumber);

	/**
	 * return true if this type is a real atom (isHeavy or H)
	 *
	 * @param type
	 * @return
	 */
	static const bool isRealAtom(const AtomTypeId & type);

	/**
	 * True if this atom type is a metal
	 *
	 * @return
	 */
	const bool isMetal() const;

	/**
	 *
	 * @return true if this is an oxygen
	 */
	static const bool isOxygenType(const AtomTypeId & type);

	/**
	 *
	 * @return true if this is a heteroatom
	 */
	static const bool isHeteroType(const AtomTypeId & type);

	/**
	 * @param type
	 *
	 * @return true if this is a halogen
	 */
	static const bool isHalogen(const AtomTypeId & type);

	/**
	 * @param type
	 *
	 * @return true if this is a carbon
	 */
	static const bool isCarbonType(const AtomTypeId & type);

	/**
	 * @param type
	 *
	 * @return true if this is a sulphur
	 */
	static const bool isSulphurType(const AtomTypeId & type);

	/**
	 * @param type
	 *
	 * @return true if this is a phosphorous
	 */
	static const bool isPhosphorousType(const AtomTypeId & type);

	/**
	 * True if this is a nitrogen atom
	 * @param type
	 * @return
	 */
	static const bool isNitrogenType(const AtomTypeId & type);

	/**
	 * Returns true if this is a heavy atom.
	 *
	 * @param type
	 * @return
	 */
	static const bool isHeavy(const AtomTypeId & type);

	/**
	 * Returns true is this is not a dummy atom or metal ion.
	 *
	 * @param type
	 * @return
	 */
	static const bool isNotDummy(const AtomTypeId & type);

	/**
	 * Tests to see if two atom types match. Checks wild card, X, Y Hev etc.
	 * Also N for N.2, N.3 etc and C, OANY, P, S.
	 *
	 * @param a2
	 * @param ignoreAromatic Set false if aromatic types are not equivalent to sp2 types (defaults to false).
	 * @return
	 */
	static const bool matchType(const AtomTypeId & type1,
			const AtomTypeId & type2, const bool ignoreAromatic = false);

	/**
	 * Tests to see if two atom types match. Checks wild card, X, Y Hev etc.
	 * Also N for N.2, N.3 etc and C, OANY, P, S. Aromatic types are not
	 * equivalent to sp2 types.
	 *
	 * @param type1
	 * @param type2
	 * @param ignoreAromatic set true for sp2 equal to aromatic, defaults to false
	 * @return
	 */
	static const bool matchElementalType(const AtomTypeId & type1,
			const AtomTypeId & type2, bool ignoreAromatic = false);

	/**
	 * Returns the elemental type. That is the non-hybridized type e.g. C for
	 * C.3
	 *
	 * @param type
	 * @return
	 */
	static const AtomTypeId getElementalType(const AtomTypeId & type);

	/**
	 * Return the SDF type string for this type.
	 *
	 * @return
	 */
	const std::string sdType() const;

	/**
	 * Returns true if this is an elemental type for which a hybridised type is
	 * available.
	 *
	 * @return
	 */
	static const bool isElementalType(const AtomTypeId & type);

	// getters and setters

	const AtomTypeId getType() const {
		return type;
	}

	const std::string& getName() const {
		return name;
	}

	const double getGaussianWidth() const {
		return gaussianWidth;
	}

	const Geometry getGeometry() const {
		return geometry;
	}

	const std::vector<int> & getNeutralBondOrders() const {
		return neutralBondOrders;
	}

	const double getRadius() const {
		return radius;
	}

	const double getWeight() const {
		return weight;
	}

	const int getAtomicNumber() const {
		return atomicNumber;
	}

	explicit AtomType(enum AtomTypeId type_, std::string name_, double radius_,
			Geometry geometry_, double weight_, std::vector<int> nbo,
			int atomicNumber_) :
			type(type_), name(name_), radius(radius_), geometry(geometry_), weight(
					weight_), neutralBondOrders(nbo), gaussianWidth(
					determineGaussianWidth(radius_)), atomicNumber(
					atomicNumber_) {
	}

	AtomType(const AtomType & rhs) = delete;
	AtomType & operator =(const AtomType & rhs) = delete;
	AtomType(AtomType && rhs) = delete;
	AtomType & operator =(AtomType && rhs) = delete;

	static const AtomType::AtomTypeId sulfonamideNitrogenType = AtomType::AtomTypeId::NPL3;

private:

	using TypesVector = std::vector< std::unique_ptr<const AtomType>>;
	const AtomTypeId type;
	const std::string name;
	const double radius;
	const Geometry geometry;
	const double weight;
	const std::vector<int> neutralBondOrders;
	const double gaussianWidth;
	const int atomicNumber;
	static constexpr const double ATOMIC_N = 2.7;

	static double determineGaussianWidth(const double radius);

	static const std::vector<AtomTypeId> nonAtomTypes;

	/**
	 * This structure holds all the static data for the class- it is created by
	 * buildTypes() and stored in atomTypes()- to prevent static initialization fiasco
	 */
	struct StaticData {

		StaticData(TypesVector types_,
				std::map<const AtomType::AtomTypeId, const AtomType *> typesLookup_,
				std::map<const std::string, const AtomType *> typesFromNameLookup_) :
				types(move(types_)), typesLookup(typesLookup_), typesFromNameLookup(
						typesFromNameLookup_) {
			;
		}

		TypesVector types;
		std::map<const AtomType::AtomTypeId, const AtomType *> typesLookup;
		std::map<const std::string, const AtomType *> typesFromNameLookup;
	};

	/**
	 * Create all atom types. Only call from atomTypes()
	 *
	 * @return
	 */
	static StaticData buildTypes();

	/**
	 * Return array of build types.
	 *
	 * @return
	 */
	static StaticData & atomTypes();

};

}

#endif /* ATOMTYPE_H_ */
