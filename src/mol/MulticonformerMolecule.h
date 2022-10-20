/*
 * MulticonformerMolecule.h
 *
 *  Created on: Aug 18, 2015
 *      Author: gjones
 */

#ifndef SRC_MOL_MULTICONFORMERMOLECULE_H_
#define SRC_MOL_MULTICONFORMERMOLECULE_H_

#include <vector>
#include "Molecule.h"
#include "mol.h"

namespace GarethMol {

using namespace std;

/**
 * Class to model a molecular conformer
 */
class Conformer {
    friend class MulticonformerMolecule;

public:
    explicit Conformer(int conformerNo_, const Molecule & mol) :
            conformerNo(conformerNo_), name(mol.getName()), coordinates(
                    mol.getCoords()) {
    }

    explicit Conformer(int conformerNo_, const Molecule & mol,
            CoordMatrix (coords)) :
            conformerNo(conformerNo_), name(mol.getName()), coordinates(coords) {
        assert(coords.cols() == static_cast<int>(mol.nAtoms()));
    }

    Conformer(const Conformer & rhs) = delete;
    Conformer & operator=(const Conformer & rhs) = delete;
    Conformer(Conformer && rhs) = delete;
    Conformer & operator=(Conformer && rhs) = delete;

    const int getConformerNo() const {
        return conformerNo;
    }

    const CoordMatrix & getCoordinates() const {
        return coordinates;
    }

    const string& getName() const {
        return name;
    }

    void setCoordinates(const CoordMatrix& coordinates) {
        this->coordinates = coordinates;
    }


    const double getEnergy() const {
        return energy;
    }

private:
    const int conformerNo;
    const string name;
    CoordMatrix coordinates;
    double energy;
};

/**
 * Represents a multiconformer molecule
 *
 */
class MulticonformerMolecule: public Molecule {

public:

    virtual ~MulticonformerMolecule() {
    }

    MulticonformerMolecule(const MulticonformerMolecule & rhs) = delete;
    MulticonformerMolecule & operator=(const MulticonformerMolecule & rhs) = delete;
    MulticonformerMolecule(const MulticonformerMolecule && rhs) = delete;
    MulticonformerMolecule & operator=(const MulticonformerMolecule && rhs) = delete;

    /**
     * Read a number of mutliconformer molecules from a file.
     *
     * @param fileName
     * @return
     */
    static std::vector<unique_ptr<MulticonformerMolecule>> readMoleculesFromFile(
            const string & fileName);

    /**
     * Read a single multiconformer molecule from a file
     *
     * @param fileName
     * @return
     */
    static unique_ptr<MulticonformerMolecule> readAndInitializeMolecule(
            const string & fileName);

    /**
     * Write multiconformer molecules to a file
     *
     * @param fileName
     * @param mols
     */
    static void writeMoleculesToFile(const string & fileName,
            const std::vector<unique_ptr<MulticonformerMolecule>> & mols);

    /**
     * Add a conformer to a molecule
     *
     * @param mol
     */
    void addConformer(const Molecule & mol) {
        conformers.emplace_back(make_unique<Conformer>(conformerNo++, mol));
    }

    const std::vector<unique_ptr<Conformer>>& getConformers() const {
        return conformers;
    }

    Conformer & getConformer(int conformerNo) {
        return *conformers.at(conformerNo);
    }

    const Conformer & getConformer(int conformerNo) const{
        return *conformers.at(conformerNo);
    }

    const size_t nConformers() const {
        return conformers.size();
    }

    int getConformerNo() const {
        return conformerNo;
    }

    /**
     * Determine Tripos Associates energies for the conformers
     */
    void calculateTaffEnergies();

private:
    std::vector<unique_ptr<Conformer>> conformers { };
    int conformerNo { 0 };

    /**
     * Constructor from base molecule.  Note at this point the input molecule is
     * invalidated as we have taken ownership of the atoms and bonds.
     *
     * @param mol
     */
    explicit MulticonformerMolecule(Molecule & mol);

    /**
     * Helper function for reading multiconformer molecules
     * @param in
     * @param fileType
     * @return
     */
    static std::vector<unique_ptr<MulticonformerMolecule>> readMolecules(
            istream & in, Molecule::FileType fileType);
};

}
/* namespace GarethMol */

#endif /* SRC_MOL_MULTICONFORMERMOLECULE_H_ */
