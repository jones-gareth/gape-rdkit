/*
 * TaffOperations.cxx
 *
 * An application for performing operations using the Tripos Associates
 * Forcefield
 *
 *  Created on: Mar 18, 2016
 *      Author: gjones
 *
 *
 */

#include "Taff.h"
#include "Molecule.h"
#include "../util/Reporter.h"

#include <boost/format.hpp>

using namespace std;
using namespace GarethMol;
using namespace GarethUtil;
using namespace GarethFF;

using MoleculePtr = unique_ptr<Molecule>;

enum class Operation {
    ENERGY, ATOM_MINIMIZE, MOL_MINIMIZE
};

/**
 * Measures energy of a series of molecules.
 *
 * @param molecules
 */
void measureEnergy(std::vector<MoleculePtr> & molecules) {
    for (auto & molecule : molecules) {

        Taff taff(*molecule);
        auto energy = taff.energy(molecule->getCoords());

        cout << "Molecule " << molecule->getName() << " Energy " << energy
                << endl;
        cout << "eVdw " << taff.getEVdw() << endl;
        cout << "eBond " << taff.getEBond() << endl;
        cout << "eAngle " << taff.getEAngle() << endl;
        cout << "eOop " << taff.getEOop() << endl;
        cout << "eTorsion " << taff.getETorsion() << endl;

    }
}

/**
 * Minimizes the energy of a series of molecules by performing atom-by-atom simplex.
 *
 * @param molecules
 * @param outFile
 */
void minimizeAtoms(std::vector<MoleculePtr> & molecules, const string outFile,
        int nCycles, double vdwDistanceCutoff, double step) {
    for (auto & molecule : molecules) {
        Taff taff(*molecule);
        if (nCycles == 0) {
            nCycles = 100;
        }
        taff.minimizeAtoms(molecule->getCoords(), nCycles, vdwDistanceCutoff, step);
    }

    Molecule::writeMoleculesToFile(outFile, molecules);
}

/**
 *
 * Minimizes the energy of a series of molecules by performing whole molecule simplex.
 *
 * @param molecules
 * @param outFile
 */
void minimizeMolecules(std::vector<MoleculePtr> & molecules,
        const string outFile, int nCycles, double step) {
    for (auto & molecule : molecules) {
        Taff taff(*molecule);
        if (nCycles == 0) {
            nCycles = 5000;
        }
        taff.minimizeMolecule(molecule->getCoords(), nCycles, step);
    }

    Molecule::writeMoleculesToFile(outFile, molecules);
}

int main(int argc, char** argv) {
    Reporter::setMinReportingLevel(Reporter::DEBUG);

    auto usage =
            (boost::format(
                    "usage: %s -e|-a|-m [-n <n_cycles>] [-c <vdw_distance_cutoff>] [-s <step>] <molecule_file_in> [<molecule_file_out>]\n")
                    % argv[0]).str();

    int c;
    auto operation = Operation::ENERGY;
    int nFiles = 1;
    int nCycles = 0;
    double vdwDistanceCutoff = 0;
    double step = 0.25;

    while ((c = getopt(argc, argv, "eamn:c:s:")) != -1) {
        switch (c) {
        case 'e':
            operation = Operation::ENERGY;
            cout << "measuring energy" << endl;
            break;
        case 'a':
            operation = Operation::ATOM_MINIMIZE;
            cout << "performing atom by atom simplex minimization" << endl;
            nFiles = 2;
            break;
        case 'm':
            operation = Operation::MOL_MINIMIZE;
            cout << "performing whole molecule simplex minimization" << endl;
            nFiles = 2;
            break;
        case 'n':
            nCycles = atoi(optarg);
            cout << "number of cycles " << nCycles << endl;
            break;
        case 'c':
            vdwDistanceCutoff = atof(optarg);
            cout << "VDW distance cutoff for atom simplex " << vdwDistanceCutoff
                    << endl;
        case 's':
            step = atof(optarg);
            cout << "Step size " << step << endl;
        }
    }

    if (argc - optind != nFiles) {
        cout << usage;
        return EXIT_FAILURE;
    }

    string molFile(argv[optind++]);
    auto molecules = Molecule::readMoleculesFromFile(molFile);
    for (auto & molecule : molecules) {
        molecule->initialize();
    }

    string outFile;
    switch (operation) {
    case Operation::ENERGY:
        measureEnergy(molecules);
        break;
    case Operation::ATOM_MINIMIZE:
        outFile = string(argv[optind]);
        minimizeAtoms(molecules, outFile, nCycles, vdwDistanceCutoff, step);
        break;
    case Operation::MOL_MINIMIZE:
        outFile = string(argv[optind]);
        minimizeMolecules(molecules, outFile, nCycles, step);
        break;
    default:
        break;
    }
}
