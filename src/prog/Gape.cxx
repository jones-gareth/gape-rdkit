//
// Created by gareth on 10/19/22.
//

#include "gape/SuperpositionMolecule.h"
#include "gape/Gape.h"
#include "Reporter.h"

#include <cassert>
#include <boost/format.hpp>
#include <GraphMol/FileParsers/MolSupplier.h>

using namespace std;
using namespace GarethUtil;

/**
 * GAPE overlap algorithm
 *
 * @param argc
 * @param argv
 * @return
 */
int main(int argc, char *argv[]) {
    assert(argc == 2);


    Reporter::setMinReportingLevel(Reporter::DEBUG);


    auto usage =
            (boost::format("usage: %s [-f] [-n <atom_no>] <smarts> <target>\n")
             % argv[0]).str();

    int c;
    bool queryFile = false;
    int atomNo = -1;

    while ((c = getopt(argc, argv, "fn:")) != -1) {
        switch (c) {
            case 'f':
                queryFile = true;
                cout << "processing file" << endl;
                break;
            case 'n':
                atomNo = atoi(optarg);
                cout << "matching smiles to atom number " << atomNo << endl;
                atomNo--;
                break;
            case '?':
                cout << usage;
                return EXIT_SUCCESS;
            default:
                cout << "Unknown option" << endl;
                cout << usage;
                return EXIT_FAILURE;
        }
    }

    if (argc - optind != 1) {
        cout << usage;
        return EXIT_FAILURE;
    }

    string inputFile(argv[optind]);
    RDKit::SmilesMolSupplier smilesMolSupplier(inputFile);
    Gape::Gape gape;
    std::vector<std::shared_ptr<Gape::SuperpositionMolecule>> molecules;
    
    while (!smilesMolSupplier.atEnd()) {
        auto mol = smilesMolSupplier.next();
        mol->setProp("_Name", "Ligand");
        auto superpositionMolecule = std::make_shared<Gape::SuperpositionMolecule>(*mol, gape.gapeSettings);
        delete mol;
        molecules.push_back(superpositionMolecule);
    }
    
    std::ofstream preparedFile("preparedMols.sdf");
    for (auto &superpositionMolecule : molecules) {
        preparedFile << superpositionMolecule->ToMolBlock() << "$$$$" << std::endl;
    }
    preparedFile.close();

}