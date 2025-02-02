//
// Created by gareth on 10/19/22.
//

#include "gape/SuperpositionMolecule.h"
#include "gape/GapeSettings.h"
#include "gape/SuperpositionGa.h"
#include "util/Reporter.h"

#include <cassert>
#include <boost/format.hpp>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/FileParsers/MolWriters.h>

#include <boost/program_options.hpp>
#include <util/GzipWriter.h>

using namespace Gape;
using namespace boost::program_options;
namespace options = boost::program_options;
using namespace RDKit;

/**
 * GAPE overlap algorithm
 *
 * @param argc
 * @param argv
 * @return
 */
int main(int argc, char *argv[]) {
    RDLog::InitLogs();
    boost::logging::disable_logs("rdApp.debug");
    Reporter::setMinReportingLevel(Reporter::DEBUG);

    std::string inputFile;
    std::string reportingLevel;
    options_description desc("Allowed options");
    desc.add_options()
            ("help", "Help message")
            ("inputFile", options::value<std::string>(&inputFile)->default_value("../../../resources/5ht3.smi"), "input structures")
            ("configFile", options::value<std::string>(), "Optional JSON configuration file")
            ("reportingLevel", options::value<std::string>(&reportingLevel)->default_value("DEBUG"),
             "Reporting level [TRACE, DEBUG, DETAIL, NORMAL, INFO, WARN, FATAL]");

    options::variables_map vm;
    options::store(options::parse_command_line(argc, argv, desc), vm);

    if (vm.count("help")) {
        cerr << desc << endl;
        return 0;
    }
    options::notify(vm);

    auto usage =
            (boost::format("usage: %s [-f] [-n <atom_no>] <smarts> <target>\n")
             % argv[0]).str();

    if (vm.count("configFile")) {

    }

	if (reportingLevel == "TRACE")
	{
        Reporter::setMinReportingLevel(Reporter::TRACE);
	}
    else if (reportingLevel == "DEBUG")
    {
        Reporter::setMinReportingLevel(Reporter::DEBUG);
    }
    else if (reportingLevel == "DETAIL")
    {
        Reporter::setMinReportingLevel(Reporter::DETAIL);
    }
    else if (reportingLevel == "NORMAL")
    {
        Reporter::setMinReportingLevel(Reporter::NORMAL);
    }
    else if (reportingLevel == "INFO")
    {
        Reporter::setMinReportingLevel(Reporter::INFO);
    }
    else if (reportingLevel == "WARN")
    {
        Reporter::setMinReportingLevel(Reporter::WARN);
    }
    else if (reportingLevel == "FATAL")
    {
        Reporter::setMinReportingLevel(Reporter::FATAL);
    }

    RDKit::SmilesMolSupplier smilesMolSupplier(inputFile, " ", 0, 1, false, true);
    GapeSettings settings;
    std::vector<std::shared_ptr<Gape::SuperpositionMolecule>> molecules;

    int ligandNum = 0;
    while (!smilesMolSupplier.atEnd()) {
        auto mol = smilesMolSupplier.next();
        ligandNum++;
        if (!mol->hasProp("_Name")) {
            mol->setProp("_Name", (boost::format("Ligand %d") % ligandNum).str());
        }
        auto superpositionMolecule = std::make_shared<Gape::SuperpositionMolecule>(*mol, settings);
        delete mol;
        superpositionMolecule->solvate();
        superpositionMolecule->generate3D();
        superpositionMolecule->findFreelyRotatableBonds();
        superpositionMolecule->findPairsToCheck();
        superpositionMolecule->findCharges();
        superpositionMolecule->findDonorsAndAcceptors();
        superpositionMolecule->findFeatures();

        molecules.push_back(superpositionMolecule);
    }

    GzipWriter gzipWriter("preparedMols.sdf.gz");
    SDWriter sdWriter(&gzipWriter.getOut());
    for (auto &superpositionMolecule: molecules) {
        sdWriter.write(superpositionMolecule->getMol());
    }
    sdWriter.close();

    Superposition superposition(molecules, settings);
    SuperpositionGa ga(superposition);
    ga.run();
}