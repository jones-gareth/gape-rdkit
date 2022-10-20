/*
 * sdfParser.cpp
 *
 *  Created on: May 17, 2014
 *      Author: Gareth Jones
 */

#include <cassert>
#include <iostream>
#include <string>
#include <vector>

#include "Reporter.h"
#include "RocsSdfParser.h"

using namespace std;
using namespace GarethUtil;
using namespace GarethMol;
using namespace Difgape;

/**
 * Test routine for reading molecules
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

	return 1;
}

