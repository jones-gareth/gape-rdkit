/*
 * MulticonformerMolecule.cpp
 *
 *  Created on: Aug 18, 2015
 *      Author: gjones
 */

#include "MulticonformerMolecule.h"

#include <iostream>
#include <memory>

#include "../util/GzipReader.h"
#include "../util/GzipWriter.h"
#include "../util/FileSystemUtil.h"
#include "Molecule.h"
#include "../util/Reporter.h"
#include "../ff/Taff.h"

namespace GarethMol {

using namespace std;
using FileType = Molecule::FileType;

	 MulticonformerMolecule::MulticonformerMolecule(Molecule & mol) :
			Molecule(mol.name, mol.atoms, mol.bonds) {
		name = mol.getName();
		sdfChiral = mol.isSdfChiral();
		conformers.reserve(250);
		// note; at this point the first conformer and the base molecule will share co-ordinates
		coords = mol.getCoords();
		addConformer(mol);
		hasCoordinates = mol.hasCoordinates;
	}

/**
 * Helper function for reading multiconformer molecules
 * @param in
 * @param fileType
 * @return
 */
vector<unique_ptr<MulticonformerMolecule>> MulticonformerMolecule::readMolecules(istream & in,
		Molecule::FileType fileType) {
	std::vector<unique_ptr<MulticonformerMolecule>> mols;

	MulticonformerMolecule * currentMultiMolecule = nullptr;
	while (true) {
		Molecule molecule;
		if (!molecule.loadMolecule(fileType, in)) {
			break;
		}

		if (currentMultiMolecule == nullptr
				|| !molecule.same2dMolecule(*currentMultiMolecule)) {
			if (currentMultiMolecule != nullptr) {
				REPORT(Reporter::NORMAL) << "loaded multiconformer molecule "
						<< currentMultiMolecule->getName() << " "
						<< to_string(currentMultiMolecule->nConformers()) + " conformers";
			}
			currentMultiMolecule = new MulticonformerMolecule(molecule);
			unique_ptr<MulticonformerMolecule> currentMultiMoleculeUP(currentMultiMolecule);
			REPORT(Reporter::DEBUG)
					<< "loaded first conformation of multiconformer molecule "
					<< currentMultiMolecule->getName();
			mols.push_back(move(currentMultiMoleculeUP));
		} else {
			currentMultiMolecule->addConformer(molecule);
			REPORT(Reporter::DEBUG) << "added conformation "
					<< currentMultiMolecule->getConformerNo()
					<< " to multiconformer molecule "
					<< currentMultiMolecule->getName();
		}
	}

	if (currentMultiMolecule != nullptr) {
		REPORT(Reporter::NORMAL) << "loaded multiconformer molecule "
				<< currentMultiMolecule->getName() << " "
				<< to_string(currentMultiMolecule->nConformers()) + " conformers";
	}
	return mols;
}

std::vector<unique_ptr<MulticonformerMolecule>> MulticonformerMolecule::readMoleculesFromFile(
		const string & fileName) {

	auto fileType = Molecule::fileNameToType(fileName);
	auto func = [fileType] (istream & in) {return readMolecules(in,fileType);};
	return readObjectFromFile<std::vector<unique_ptr<MulticonformerMolecule>>>(fileName, func);
}

unique_ptr<MulticonformerMolecule> MulticonformerMolecule::readAndInitializeMolecule(
		const string & fileName) {
	auto molecules = readMoleculesFromFile(fileName);
	assert(molecules.size() == 1ul);
	auto molecule = move(molecules.at(0));
	molecule->initialize();
	return molecule;
}

/**
 * Helper function for writing multiconformer molecules
 * @param out
 * @param mols
 * @param fileType
 */
static void writeMolecules(ostream & out,
		const std::vector<unique_ptr<MulticonformerMolecule>> & mols,
		const Molecule::FileType & fileType) {
	for (auto &  molecule : mols) {
		for (auto & conformer : molecule->getConformers()) {
			auto & otherCoords = conformer->getCoordinates();
			molecule->writeToFile(fileType, out, otherCoords, false, "");
		}
	}
}

void MulticonformerMolecule::writeMoleculesToFile(const string & fileName,
		const std::vector<unique_ptr<MulticonformerMolecule>> & mols) {

	auto fileType = Molecule::fileNameToType(fileName);
	auto func =
			[& mols, fileType](ostream & out) {writeMolecules(out, mols, fileType);};
	writeObjectToFile(fileName, func);

}

	/**
	 * Calculate TAFF forcefield energies for each conformer in the molecule
	 */
	void MulticonformerMolecule::calculateTaffEnergies() {
        GarethFF::Taff taff(*this);
        for (auto & conformer: conformers) {
            auto energy = taff.energy(conformer->getCoordinates());
            conformer->energy = energy;
        }
	}
} /* namespace GarethMol */
