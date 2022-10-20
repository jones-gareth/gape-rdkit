/*
 * MatchSmartsToStructure.cxx
 *
 *  Created on: Oct 19, 2015
 *      Author: gjones
 */

#include <iostream>
#include <string>

#include <boost/format.hpp>
#include "../mol/SmilesParser.h"
#include "../mol/SmartsParser.h"
#include "../mol/SubstructureSearch.h"
#include "../util/Reporter.h"

using namespace std;
using namespace GarethMol;
using namespace GarethUtil;

/**
 * Application that matches a smarts pattern to a structure.
 *
 * @param argc
 * @param argv
 * @return
 */
int main(int argc, char* argv[]) {

	Reporter::setMinReportingLevel(Reporter::DETAIL);

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

	if (argc - optind != 2) {
		cout << usage;
		return EXIT_FAILURE;
	}

	string smarts(argv[optind]);
	string target(argv[optind + 1]);

	SmartsParser smartsParser;
	SmilesParser smilesParser;
	unique_ptr<Molecule> targetMolecule = nullptr;
	if (queryFile) {
		// target is a file name
		targetMolecule = Molecule::readAndInitializeMolecule(target);
	} else {
		// target is a smiles string
		targetMolecule = smilesParser.parseSmiles(target);
		ofstream out("test_smiles.sdf");
		targetMolecule->writeSdfFile(out, "molecule from smiles " + target);
		targetMolecule->initialize();
	}

	auto queryMolecule = smartsParser.parseSmarts(smarts);

	if (atomNo < 0) {
		SubstructureSearch search(*queryMolecule, *targetMolecule);

		auto callbackFunction =
				[& targetMolecule] (const vector<size_t> & queryIdsToTargetIds) {
					cout << "Found isomorphism"<<endl;
					for (size_t queryId = 0; queryId < queryIdsToTargetIds.size(); queryId++) {
						auto targetId = queryIdsToTargetIds.at(queryId);
						auto & targetAtom = targetMolecule->getAtom(targetId);
						cout << "Query atom " << queryId << " maps to target atom "
						<< targetAtom.info() << endl;
					}
				};
		search.setCallbackFunction(callbackFunction);

		auto nIsomorphisms = search.compare();
		cout << "Found " << nIsomorphisms << " isomorphisms";
	} else {
		const auto & atom = targetMolecule->getAtom(atomNo);
		bool match = atom.matchAtom(smarts, *targetMolecule);
		if (match) {
			cout << "Atom " << to_string(atomNo + 1)
					<< " matches first atom in smarts pattern";
		} else {
			cout << "Atom " << to_string(atomNo + 1)
					<< " does not match first atom in smarts pattern";
		}
	}

	return EXIT_SUCCESS;
}
