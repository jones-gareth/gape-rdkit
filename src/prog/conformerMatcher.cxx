/*
 * conformerMatcher.cpp
 *
 *  Created on: Jun 12, 2014
 *      Author: Gareth Jones
 */

#include <cassert>
#include <iostream>
#include <string>
#include <vector>

#include "Reporter.h"
#include "RocsSdfParser.h"
#include "ConformerMatcher.h"

using namespace std;
using namespace GarethUtil;
using namespace Difgape;

/**
 * Difgape overlap algorithm
 *
 * @param argc
 * @param argv
 * @return
 */
int main(int argc, char* argv[]) {
	assert(argc == 2);
	string directory {argv[1]};

	Reporter::setMinReportingLevel(Reporter::INFO);

	RocsSdfParser parser {};
	parser.readDirectory(directory);

	ConformerMatcher conformerMatcher {parser};
	conformerMatcher.run();

	return 1;
}


