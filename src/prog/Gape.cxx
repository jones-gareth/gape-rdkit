//
// Created by gareth on 10/19/22.
//

#include "gape/SuperpositionMolecule.h"
#include "gape/GapeSettings.h"
#include "gape/SuperpositionGa.h"
#include "util/Reporter.h"

#include <cassert>
#include <filesystem>
#include <boost/format.hpp>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/FileParsers/MolWriters.h>

#include <boost/program_options.hpp>
#include <util/GzipWriter.h>

using namespace Gape;
using namespace boost::program_options;
namespace options = boost::program_options;
using namespace RDKit;

namespace DETAIL {
    void cleanupFiles() {
        for (const auto& entry: std::filesystem::directory_iterator(".")) {
            const auto& filename = entry.path().filename().string();
            std::string cmp = "GA_rank_";
            std::string cmp2 = "GA_solution_";
            if (filename.compare(0, cmp.length(), cmp) == 0) {
                std::filesystem::remove(entry.path());
            }
            else if (filename.compare(0, cmp2.length(), cmp2) == 0) {
                std::filesystem::remove(entry.path());
            }
        }
    }
}

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
    std::string configFile;
    // options_description desc("Allowed options");
    options_description desc;
    desc.add_options()
            ("help", "Help message")
            ("inputFile", options::value<std::string>(&inputFile)->default_value("../../../resources/5ht3.smi"), "input structures")
            ("configFile", options::value<std::string>(&configFile), "Optional JSON configuration file")
            ("reportingLevel", options::value<std::string>(&reportingLevel)->default_value("DEBUG"),
             "Reporting level [TRACE, DEBUG, DETAIL, NORMAL, INFO, WARN, FATAL]");

    options::variables_map vm;
    options::store(options::parse_command_line(argc, argv, desc), vm);

    if (vm.count("help")) {
        cerr << desc << endl;
        return 0;
    }
    options::notify(vm);

    DETAIL::cleanupFiles();

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

    GapeSettings settings(configFile);
    if (false) {
        // show that RDKIt is not responsible for the parallel slowness
        SuperpositionGa::testBatchRun(inputFile, settings);
    } else {
        RDKit::SmilesMolSupplier smilesMolSupplier(inputFile, " ", 0, 1, false, true);
        auto molecules = SuperpositionMolecule::loadMolecules(smilesMolSupplier, settings);

        {
            GzipWriter gzipWriter("preparedMols.sdf.gz");
            SDWriter sdWriter(&gzipWriter.getOut());
            for (auto &superpositionMolecule: molecules) {
                sdWriter.write(superpositionMolecule->getMol());
            }
            sdWriter.close();
        }


        Superposition superposition(molecules, settings);
        SuperpositionGa::batchRun(superposition);
    }
}