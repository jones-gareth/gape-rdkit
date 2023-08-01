/* 
 * File:   testChrom.cpp
 * Author: Gareth Jones
 *
 * Created on April 12, 2013, 8:52 AM
 */

#include <cstdlib>
#include <cstdint>
#include <ctime>
#include <iostream>
#include <vector>
#include <map>
#include <stdexcept>
#include <boost/format.hpp>
#include <cmath>

#include "util/Reporter.h"
#include "util/RandomUtil.h"
#include "ga/BinaryStringChromosomePolicy.h"
#include "ga/IntegerStringChromosomePolicy.h"
#include "ga/StringChromosome.h"
#include "ga/LinkedPopLinearSel.h"

using namespace std;
using namespace Gape;

void testBinaryChrom(RandomUtil & rng);
void testIntegerChrom(RandomUtil & rng);

/**
 * Test routine for chromosomes
 * 
 * @param argc
 * @param argv
 * @return
 */
int main(int argc, char** argv) {

	uint32_t seed = time(nullptr);
	auto & rng = RandomUtil::getInstance();
	rng.seed(seed);

	cout << "showing binary ops " << endl;
	testBinaryChrom(rng);
	cout << endl << "showing integer ops " << endl;
	testIntegerChrom(rng);

	cout << boost::format("writing a double %|10.2f|\n") % 1.23456789;

	std::multimap<double, std::string> test;

	test.insert(std::pair<double, std::string>(134, "first"));
	test.insert(std::pair<double, std::string>(1.234, "hello"));
	test.insert(std::pair<double, std::string>(1.234, "hello 2"));
	test.insert(std::pair<double, std::string>(-134, "last"));

	cout << "searching list" << endl;
	auto iterators = test.equal_range(-134);
	for (auto iterator = iterators.first; iterator != iterators.second; ++iterator) {
		cout << iterator->second << endl;
	}
}

void testBinaryChrom(RandomUtil & rng) {

	const int length = 20;

	BinaryStringChromosomePolicy chromPolicy(rng);
	BinaryStringChromosome p1(length, rng, chromPolicy);
	p1.initialize();
	BinaryStringChromosome p2(length, rng, chromPolicy);
	p2.initialize();
	BinaryStringChromosome c1(length, rng, chromPolicy);
	BinaryStringChromosome c2(length, rng, chromPolicy);
	p1.onePointCrossover(p2, c1, c2);

	cout << "Crossover" << endl;
	cout << "Parent 1| " << p1.geneInfo() << endl;
	cout << "Parent 2| " << p2.geneInfo() << endl;
	cout << "Child  1| " << c1.geneInfo() << endl;
	cout << "Child  2| " << c2.geneInfo() << endl;
	cout << " test decode " << c2.decodeByte(1) << endl;

	BinaryStringChromosome mp(length, rng, chromPolicy);
	mp.initialize();
	BinaryStringChromosome mc(length, rng, chromPolicy);
	mc.copyGene(mp);
	cout << "Mutation" << endl;
	cout << "Parent  | " << mp.geneInfo() << endl;
	mc.mutate();
	cout << "Child   | " << mc.geneInfo() << endl;

}

void testIntegerChrom(RandomUtil & rng) {

	const int length = 20;
	IntegerStringChromosomePolicy chromPolicy(rng, 20);
	chromPolicy.setMax(10);
	chromPolicy.setAllowNulls(false);
	IntegerStringChromosome p1(length, rng, chromPolicy);
	p1.initialize();
	IntegerStringChromosome p2(length, rng, chromPolicy);
	p2.initialize();
	IntegerStringChromosome c1(length, rng, chromPolicy);
	IntegerStringChromosome c2(length, rng, chromPolicy);
	p1.fullMixingAndCrossover(p2, c1, c2);

	cout << "Crossover" << endl;
	cout << "Parent 1| " << p1.geneInfo() << endl;
	cout << "Parent 2| " << p2.geneInfo() << endl;
	cout << "Child  1| " << c1.geneInfo() << endl;
	cout << "Child  2| " << c2.geneInfo() << endl;

	IntegerStringChromosome mp(length, rng, chromPolicy);
	mp.initialize();
	IntegerStringChromosome mc(length, rng, chromPolicy);
	mc.copyGene(mp);
	cout << "Mutation" << endl;
	cout << "Parent  | " << mp.geneInfo() << endl;
	mc.mutate();
	cout << "Child   | " << mc.geneInfo() << endl;

}
