/*
 * testGa.cpp
 *
 * A test program for a binary GA
 *
 *  Created on: Apr 26, 2013
 *      Author: gjones
 */

#include "Reporter.h"
#include "BinaryTestGa.h"

using namespace GapeGa;

int main(int argc, char* argv[]) {
    Reporter::setMinReportingLevel(Reporter::TRACE);
    BinaryTestGa ga;
    ga.run();
    return 1;
}

