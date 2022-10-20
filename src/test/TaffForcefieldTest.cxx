/*
 * TaffForcefieldTest
 *
 *	Tests the Tripos force field.
 *
 *	Clark, M. . Cramer., R. D. .. Van Opdenbosch,N. Validation of the General Purpose
 *	Tripos 5.2 Force Field. J.Comp.Chem. 10, 982–1012 (1989).
 *
 */

#include <gtest/gtest.h>
#include <boost/range/irange.hpp>
// #include "TestUtil.h"
#include "../ff/Taff.h"
#include "../util/Reporter.h"

namespace {

using namespace std;
using namespace GarethUtil;
// using namespace GarethTest;
using namespace GarethFF;

using namespace ::testing;

class TaffForcefieldTest: public Test {

protected:

    static void SetUpTestCase() {
        Reporter::setMinReportingLevel(Reporter::DETAIL);
        //Reporter::setFileStream(out);
    }

    static void TearDownTestCase() {
        out.close();
    }

    static ofstream out;
};

ofstream TaffForcefieldTest::out { "TaffForcefieldTest.log" };

/**
 * Reads all ligands from the Astex diverse set and checks their energy.
 *
 * The reference energies come from my Java code.  So this checks that the programs
 *  are consistent. I don't have access to SYBL so I don't know how to produce exact
 *  reference energies.
 */
TEST_F(TaffForcefieldTest, TestAllLigands) {
    std::vector<double> testEnergies = { 57.86884577169073, 79.40985071692187,
            164.6276675197239, 61.75851212080544, 67.08984840414007,
            157.13663641577122, 142.4541330657746, 186.04370394883927,
            159.32312982727646, 75.50253954945673, 83.44189716601906,
            124.78565541526909, 75.93659766620496, 106.72692380164584,
            61.73177219740269, 86.51311140005653, 68.3316505224013,
            284.8123226851252, 110.50094541626655, 269.55020065559194,
            83.36765577110329, 255.2183278064647, 41.132985928491266,
            116.79268216208933, 75.99564339923431, 106.23896691800651,
            134.97022979644086, 168.8802443591628, 70.7265488211068,
            121.21003345056552, 110.29251406940742, 133.01682667557682,
            520.1017643263164, 38.78542379222401, 156.47421045584008,
            110.28438858080982, 103.85288966029945, 154.09002732850436,
            37.9960651705849, 129.2833874818305, 497.48600803349655,
            158.84306875242797, 166.99484065958666, 37.941057725927145,
            65.79944405072129, 65.97351220789055, 158.48904648258255,
            95.33815365902771, 505.08700587903473, 95.23734631659939,
            90.47703772327316, 128.84126170374734, 100.03027859472677,
            153.64571057604928, 153.01995026945497, 84.48561412556396,
            217.01839365358083, 73.71031657674006, 61.56004387885703,
            62.02911407645805, 145.5530735864179, 249.32896673047526,
            161.57075380320572, 95.70287543345337, 55.727563427544695,
            185.12610270136486, 186.74070927860598, 158.69906300139016,
            112.46734284303461, 39.12331403237114, 108.75130064044424,
            25.295516988612885, 60.35538291520852, 319.5857324782783,
            113.6891411987643, 323.86577109794297, 253.341950440147,
            112.71932625482579, 126.41610034284493, 131.48579761583113,
            104.6272273421172, 101.36819957467787, 102.92636808962354,
            183.34477372032933, 91.5594765603172 };

    auto molecules = Molecule::readMoleculesFromFile(
            "../resources/all_ligands_charged.mol");
    ASSERT_EQ(molecules.size(), testEnergies.size());

    for (auto i : boost::irange(0ul, molecules.size())) {
        auto &molecule = molecules.at(i);
        REPORT(Reporter::DETAIL) << "Testing ligand " << molecule->getName();
        molecule->initialize();
        auto expectedEnergy = testEnergies.at(i);
        Taff taff(*molecule);
        auto energy = taff.energy(molecule->getCoords());
        REPORT(Reporter::DETAIL) << "Ligand no " << to_string(i + 1) << " "
                << molecule->getName() << " expected energy " << expectedEnergy
                << " calculated energy " << energy;
        ASSERT_TRUE(equals(expectedEnergy, energy, 0.0001));
    }

}

}
