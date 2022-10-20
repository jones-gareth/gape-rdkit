/*
 * loadMol.cpp
 *
 *	Test file for reading in smiles
 *
 *  Created on: Apr 29, 2014
 *      Author: gjones
 */

#include <gtest/gtest.h>
#include <stddef.h>
#include <cstdio>
#include <iostream>
#include <memory>
#include <string>
#include <vector>
#include <fstream>

#include "../mol/Molecule.h"
#include "../mol/MolComparer.h"
#include "../mol/SmartsParser.h"
#include "../util/Reporter.h"

namespace {

using namespace std;
using namespace GarethUtil;
using namespace GarethMol;
using namespace ::testing;

/**
 * Test class for reading samrts patterns
 */
class SmartsTest: public Test {

protected:
	SmartsTest() {
		Reporter::setMinReportingLevel(Reporter::DETAIL);
	}
};

TEST_F(SmartsTest, testEquiv1) {
	auto smiles1 = "CC1=CC(Br)CCC1";
	auto smiles2 = "CC1=CC(CCC1)Br";
	SmartsParser parser;
	parser.parseSmarts(smiles1);
	parser.parseSmarts(smiles2);
}

//void test2 () {
TEST_F(SmartsTest, testCanLoad) {
	ifstream in("../resources/smarts_library.ok.test");
	ASSERT_TRUE(in.is_open());
	SmartsParser parser;
	auto cnt = 0;
	while (!in.eof()) {
		string line;
		getline(in, line);
		removeTrailingLF(line);
		trim(line);
		if (line == "") {
			continue;
		}
		if (startsWith(line, "#")) {
			continue;
		}
		cnt++;
		REPORT(Reporter::DETAIL) << "Testing parse of smarts no " << cnt
				<< ": " << line;
		parser.parseSmarts(line);
	}
}

}
