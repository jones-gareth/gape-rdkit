/*
 * testGa.cpp
 *
 * A test program for a binary GA
 *
 *  Created on: Apr 26, 2013
 *      Author: Gareth Jones
 */

#include "../util/Reporter.h"
#include "../ga/BinaryTestGa.h"

using namespace Gape;

int main(int argc, char* argv[]) {
    Reporter::setMinReportingLevel(Reporter::TRACE);
    BinaryTestGa ga;
    ga.run();
    return 1;
}

