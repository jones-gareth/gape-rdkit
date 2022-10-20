/*
 * Molecule.h

 *
 *  Created on: Apr 28, 2014
 *      Author: gjones
 */

#ifndef MOLECULE_H_
#define MOLECULE_H_

#include <stddef.h>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>
#include <boost/optional.hpp>
#include <Eigen/Dense>

#include "../util/Util.h"
#include "../util/CoordOps.h"
#include "Atom.h"
#include "Bond.h"
#include "mol.h"
#include "Ring.h"
#include "../util/Array2D.h"

namespace GarethMol {

using namespace std;
using namespace GarethUtil;
using namespace Eigen;

class Torsion {
public:
	Torsion(const int no, const int a1, const int a2, const int a3, const int a4, const int b) :
			torsionNo(no), atom1(a1), atom2(a2), atom3(a3), atom4(a4), bond23(b) {
	}

	virtual ~Torsion() {
	}

	Torsion(const Torsion & rhs) = delete;
	Torsion & operator =(const Torsion & rhs) = delete;
	Torsion(Torsion && rhs) = delete;
	Torsion & operator =(Torsion && rhs) = delete;

	const int getAtom1() const {
		return atom1;
	}

	const int getAtom2() const {
		return atom2;
	}

	const int getAtom3() const {
		return atom3;
	}

	const int getAtom4() const {
		return atom4;
	}

    const int getBond23() const {
        return bond23;
    }

	const int getTorsionNo() const {
		return torsionNo;
	}

private:
	const int torsionNo, atom1, atom2, atom3, atom4, bond23;
};

class Angle {
public:
	Angle(const int no, const int a1, const int a2, const int a3) :
			angleNo(no), atom1(a1), atom2(a2), atom3(a3) {
	}

	virtual ~Angle() {
	}

	Angle(const Angle & rhs) = delete;
	Angle & operator =(const Angle & rhs) = delete;
	Angle(Angle && rhs) = delete;
	Angle & operator =(Angle && rhs) = delete;

	const int getAtom1() const {
		return atom1;
	}

	const int getAtom2() const {
		return atom2;
	}

	const int getAtom3() const {
		return atom3;
	}

	const int getAngleNo() const {
		return angleNo;
	}

private:
	const int angleNo, atom1, atom2, atom3;
};
/**
 * Class for a molecule
 */
class Molecule {
	friend class AtomAddition;
	friend class MulticonformerMolecule;

public:
	using MoleculePtr = unique_ptr<Molecule>;

	enum class FileType {
		SDF, MOL2
	};

	explicit Molecule() {
	}

	explicit Molecule(const string & name_,
			std::vector<unique_ptr<Atom>> & atoms_,
			std::vector<unique_ptr<Bond>> & bonds_) :
			name(name_), atoms(move(atoms_)), bonds(move(bonds_)) {
		coords = CoordMatrix::Zero(4, atoms.size());
		coords.row(3) = VectorXd::Constant(atoms.size(), 1.0);
		//cout << "Here is the coordinate matrix:\n" << coords<<endl;
		hasCoordinates = false;
	}

	virtual ~Molecule() {
	}

	Molecule(const Molecule & rhs) = delete;
	Molecule & operator =(const Molecule & rhs) = delete;
	Molecule(Molecule && rhs) = delete;
	Molecule & operator =(Molecule && rhs) = delete;

	/**
	 * load a MOL2 molecule and return it.
	 *
	 * @param in
	 * @return
	 */
	bool loadMol2Molecule(istream & in);

	/**
	 * Load an SDF molecule and return it.
	 *
	 * @param in
	 * @return
	 */
	bool loadSdfMolecule(istream & in);

	/**
	 * Load a molecule of the appropriate file type
	 *
	 * @param fileType
	 * @param in
	 * @return
	 */
	bool loadMolecule(const FileType fileType, istream &in);

	/**
	 * Write the molecule to the file using the format specified.
	 *
	 * @param fileType
	 * @param out
	 * @param includeLonePairs
	 * @param comment
	 */
	void writeToFile(const FileType & fileType, ostream & out,
			const bool & includeLonePairs = true,
			const string & comment = "") const;

	/**
	 * Write the molecule to the file using the format and coordinates specified.
	 * @param fileType
	 * @param out
	 * @param otherCoords
	 * @param includeLonePairs
	 * @param comment
	 */
	void writeToFile(const FileType & fileType, ostream & out,
			const CoordMatrix & otherCoords, const bool & includeLonePairs =
					true, const string & comment = "") const;

	/**
	 * write the molecule to the output stream in MOL2 format.
	 *
	 * @param out
	 * @param comment
	 */
	void writeMol2File(ostream & out, const bool & includeLonePairs = true,
			const string & comment = "") const;

	/**
	 * Write the molecule out with alternative coordinates
	 * @param out
	 * @param otherCoords
	 * @param includeLonePairs
	 * @param comment
	 */
	void writeMol2File(ostream & out, const CoordMatrix & otherCoords,
			const bool & includeLonePairs = true,
			const string & comment = "") const;

	/**
	 * write the molecule to the output stream in SDF format.
	 *
	 * @param out
	 * @param comment
	 */
	void writeSdfFile(ostream & out, const string & comment = "") const;

	/**
	 * Write the molecule out with alternative coordinates
	 *
	 * @param out
	 * @param otherCoords
	 * @param comment
	 */
	void writeSdfFile(ostream & out, const CoordMatrix & otherCoords,
			const string & comment = "") const;

	/**
	 * Determine the file type of a file name
	 * @param fileName
	 * @param gzipped return true if gzipped
	 * @return
	 */
	static FileType fileNameToType(const string & fileName);

	/**
	 * Returns true if this molecule has the same connection table as the other molecule.
	 * Does not account for isomorphisms (yet).
	 *
	 * @param otherMol
	 * @return
	 */
	bool same2dMolecule(const Molecule & otherMol) const;

	/**
	 * @return all atoms in the molecule
	 */
	const std::vector<unique_ptr<Atom>>& getAtoms() const {
		return atoms;
	}

	/**
	 * @return all bonds in the molecule
	 */
	const std::vector<unique_ptr<Bond>>& getBonds() const {
		return bonds;
	}

	/**
	 * @return all molecular coordinates
	 */
	virtual const CoordMatrix getCoords() const {
		return coords;
	}

    /**
     * @return all molecular coordinates (for editing)
     */
    CoordMatrix & getCoords()  {
        return coords;
    }

	const CoordVector getCoord(int atomNo) const {
		return coords.col(atomNo);
	}

	/**
	 * @param atomNo
	 * @return A particular atom
	 */
	Atom & getAtom(int atomNo) const {
		return *atoms.at(atomNo);
	}

	/**
	 * @param bondNo
	 * @return a particular bond
	 */
	Bond & getBond(int bondNo) const {
		return *bonds.at(bondNo);
	}

	/**
	 * @return number of atoms in the molecule
	 */
	size_t nAtoms() const {
		return atoms.size();
	}

	/**
	 * @return number of bonds in the molecule
	 */
	size_t nBonds() const {
		return bonds.size();
	}

	/**
	 * @return molecular name
	 */
	const string& getName() const {
		return name;
	}

	/**
	 * Set molecular name
	 *
	 * @param name
	 */
	void setName(const string& name) {
		this->name = name;
	}

	/**
	 * Add a sdf field to the molecule
	 *
	 * @param field
	 * @param values
	 */
	void addSdfField(string field, std::vector<string> values) {
		sdfFields[field] = values;
	}

	/**
	 * Get the values for an sd field
	 *
	 * @param field
	 * @return
	 */
	std::vector<string> getSdfValues(const string & field) const {
		return sdfFields.at(field);
	}

	/**
	 * Checks if an sdf field is present
	 *
	 * @param field
	 * @return
	 */
	bool hasSdfField(const string & field) const {
		return isKeyPresent(field, sdfFields);
	}

	/**
	 *
	 * @return the map of all sd fields and values
	 */
	const map<string, std::vector<string> >& getSdfFields() const {
		return sdfFields;
	}

	/**
	 * Reads the molecule from a file.
	 *
	 * @param fileName
	 */
	bool readMoleculeFromFile(const string & fileName);

	static std::vector<MoleculePtr> readMoleculesFromFile(
			const string & fileName);

	static unique_ptr<Molecule> readAndInitializeMolecule(
			const string & fileName);

	static void writeMoleculesToFile(string fileName,
			std::vector<MoleculePtr> & mols);

	bool isSdfChiral() const {
		return sdfChiral;
	}

	/**
	 * Return true if two atoms are bound- this searches the bond connection table
	 * @param atom1
	 * @param atom2
	 * @return
	 */
	boost::optional<const Bond &> isBonded(const Atom & atom1,
			const Atom & atom2) const;

	/**
	 * Return true if two atoms are bound- this searches the bond connection table
	 */
	boost::optional<const Bond &> isBonded(const int atom1No,
			const int atom2No) const;

	/**
	 * Create atom neighbourhood information
	 */
	void createAtomNeighbourhood();

	/**
	 * Creates atom neighbourhood and performs SSSR determination
	 */
	void initialize();

	/**
	 * Retrieve a bond from the bond lookup table.  Must call createBondTable() first.
	 *
	 * @param atom1
	 * @param atom2
	 * @return
	 */
	Bond * getBond(int atom1, int atom2) const;

	/**
	 * Return true if these atoms are seperated by two bonds (or form the ends
	 * of an angle). Must call createBondTable() first.
	 *
	 * @param atom1
	 * @param atom2
	 * @return
	 */
	bool is13Bonded(int atom1, int atom2) const{
		return bonds13->get(atom1, atom2);
	}

	/**
	 * Return true if these atoms are seperated by three bonds (or form the
	 * ends of a torsion).  Must call createBondTable() first.
	 *
	 * @param atom1
	 * @param atom2
	 * @return
	 */
	bool is14Bonded(int atom1, int atom2) const {
		return bonds14->get(atom1, atom2);
	}

	/**
	 * Return true if the two atoms are bonded using the lookup table.  Must call createBondTable() first.
	 *
	 * @param atom1No
	 * @param atom2No
	 * @return
	 */
	bool is12Bonded(const int atom1No, const int atom2No) const {
		return getBond(atom1No, atom2No) != nullptr;
	}

	void clearRings() {
		rings.clear();
	}

	void addRing(unique_ptr<Ring> & ring) {
		rings.push_back(move(ring));
	}

	const std::vector<unique_ptr<Ring> >& getRings() const {
		return rings;
	}

	/**
	 * Assign atom and bond types- returns number of atom types changed.
	 *
	 * The molcule needs a bond table, atom neighbourhoods, SSSR, before we
	 * call this.
	 *
	 * @return
	 */
	int assignAtomTypes();

	/**
	 * Centers the coordinate block
	 */
	void centerCoords();

	/**
	 * Returns the centroid of the current coordinates
	 *
	 * @return
	 */
	CoordVector centroid() const;

	/**
	 * Returns the number of real atoms (not LP etc)
	 *
	 * @return
	 */
	int nRealAtoms() const;

    /**
     * Returns the number of heavy atoms
     *
     * @return
     */
    int nHeavyAtoms() const;

	/**
	 * Returns the number of real bonds (not to LP etc)
	 *
	 * @return
	 */
	int nRealBonds() const;

	/**
	 * Creates the bond table lookup only if it's absent
	 */
	void createBondTableIfAbsent() {
		if (bondTable == nullptr) {
			createBondTable();
		}
	}

	bool isHasCoordinates() const {
		return hasCoordinates;
	}

	void setHasCoordinates(bool hasCoordinates = false) {
		this->hasCoordinates = hasCoordinates;
	}

    const std::vector<unique_ptr<Angle> >& getAngles() const {
        return angles;
    }

    const std::vector<unique_ptr<Torsion> >& getTorsions() const {
        return torsions;
    }

protected:
	string name;
	// bonds, atoms, coords etc
	std::vector<unique_ptr<Atom>> atoms { };
	std::vector<unique_ptr<Bond>> bonds { };
	std::vector<unique_ptr<Ring>> rings { };
	CoordMatrix coords;
	unique_ptr<Array2D<Bond *>> bondTable = nullptr;
	unique_ptr<Array2D<bool>> bonds13 = nullptr;
	unique_ptr<Array2D<bool>> bonds14 = nullptr;
	bool hasCoordinates = false;
	std::vector<unique_ptr<Torsion>> torsions;
	std::vector<unique_ptr<Angle>> angles;

	/**
	 * creates a mapping for atom numbers when outputting a structure and skipping lone pairs.
	 *
	 * @param includeLonePairs
	 * @return
	 */
	const map<int, int> mapAtomsForOutput(
			const bool & includeLonePairs = false) const;

	/**
	 * Creates a bond table lookup, and angle and torsion lists
	 */
	void createBondTable();

	/**
	 * Updates atom neighbourhood information with rings.
	 *
	 */
	void updateNeighbourHoodWithRings();

	// sd flags
	bool sdfChiral = false;
	//double sdMassDifference = 0.0;
	map<string, std::vector<string>> sdfFields;

};

} /* namespace GarethMol */

#endif /* MOLECULE_H_ */
