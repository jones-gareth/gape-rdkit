/*
 * LeastSquaresTest.cpp
 *
 *  Created on: Sep 10, 2015
 *      Author: gjones
 */

#include <gtest/gtest.h>
#include "../mol/Molecule.h"
#include "../mol/MolComparer.h"
#include "../mol/MolecularRms.h"
#include "../mol/MulticonformerMolecule.h"
#include "../mol/MultiConformerRmsEvalution.h"

namespace {

using namespace std;
using namespace GarethUtil;
using namespace GarethMol;
using namespace ::testing;

class LeastSquaresTest: public Test {

protected:
	LeastSquaresTest() {
		Reporter::setMinReportingLevel(Reporter::NORMAL);
	}

};

static const double TOLERANCE = 1e-4;

// Reference values for isomorphic rms least squares between 20 conformers
// of 1h9v from java program.
static const map<pair<size_t, size_t>, double> omegaResults { { make_pair(0, 1),
		1.1472572773227063 }, { make_pair(0, 2), 1.1749003032105503 }, {
		make_pair(0, 3), 1.4797540569663412 }, { make_pair(0, 4),
		0.7385902758642684 }, { make_pair(0, 5), 1.1101221035249746 }, {
		make_pair(0, 6), 1.1560661859579748 }, { make_pair(0, 7),
		1.5040988503082136 }, { make_pair(0, 8), 1.0991767581584748 }, {
		make_pair(0, 9), 1.2083880244212633 }, { make_pair(0, 10),
		1.3887264671502313 }, { make_pair(0, 11), 1.3737233211846775 }, {
		make_pair(0, 12), 1.0918741575278224 }, { make_pair(0, 13),
		1.4449914204231706 }, { make_pair(0, 14), 1.3507289253427208 }, {
		make_pair(0, 15), 1.3579358371949803 }, { make_pair(0, 16),
		1.2973157972442098 }, { make_pair(0, 17), 1.6882912442824491 }, {
		make_pair(0, 18), 1.0007577054706476 }, { make_pair(0, 19),
		1.5215195544008577 }, { make_pair(1, 2), 1.4814157147221283 }, {
		make_pair(1, 3), 1.17450904481589 }, { make_pair(1, 4),
		1.5039765987154543 }, { make_pair(1, 5), 1.1552342052779188 }, {
		make_pair(1, 6), 1.110625596768799 }, { make_pair(1, 7),
		0.7387776641490573 }, { make_pair(1, 8), 1.2072476817410644 }, {
		make_pair(1, 9), 1.0991615226552875 }, { make_pair(1, 10),
		1.3734246459503958 }, { make_pair(1, 11), 1.3898042449175694 }, {
		make_pair(1, 12), 1.44426197552555 }, { make_pair(1, 13),
		1.0911064420808354 }, { make_pair(1, 14), 1.3577084693810861 }, {
		make_pair(1, 15), 1.3509192402986274 }, { make_pair(1, 16),
		1.6553725406017787 }, { make_pair(1, 17), 1.2756799195722286 }, {
		make_pair(1, 18), 1.4107775949289185 }, { make_pair(1, 19),
		1.0215140126370197 }, { make_pair(2, 3), 1.125680092413449 }, {
		make_pair(2, 4), 1.146519778896903 }, { make_pair(2, 5),
		1.3453876710310133 }, { make_pair(2, 6), 0.7750805742513119 }, {
		make_pair(2, 7), 1.1834216649760299 }, { make_pair(2, 8),
		1.7249786040467558 }, { make_pair(2, 9), 1.8098277176353412 }, {
		make_pair(2, 10), 0.9994054843790348 }, { make_pair(2, 11),
		1.222481493879754 }, { make_pair(2, 12), 1.7849366046704416 }, {
		make_pair(2, 13), 1.7773687652954495 }, { make_pair(2, 14),
		1.3391539128590209 }, { make_pair(2, 15), 1.0984064606372526 }, {
		make_pair(2, 16), 1.5316163026691563 }, { make_pair(2, 17),
		1.8356159316083258 }, { make_pair(2, 18), 1.4750763430120635 }, {
		make_pair(2, 19), 1.656621651557799 }, { make_pair(3, 4),
		1.1813102176767976 }, { make_pair(3, 5), 0.7741762665470789 }, {
		make_pair(3, 6), 1.343824057188411 }, { make_pair(3, 7),
		1.1467598394522702 }, { make_pair(3, 8), 1.8076690532611879 }, {
		make_pair(3, 9), 1.7237830407760466 }, { make_pair(3, 10),
		1.2207958865801494 }, { make_pair(3, 11), 0.9990026939852976 }, {
		make_pair(3, 12), 1.7749987675926964 }, { make_pair(3, 13),
		1.7839554750363191 }, { make_pair(3, 14), 1.0975287133769795 }, {
		make_pair(3, 15), 1.3372489568932133 }, { make_pair(3, 16),
		1.7561654243801397 }, { make_pair(3, 17), 1.4564300473130873 }, {
		make_pair(3, 18), 1.7986589884805377 }, { make_pair(3, 19),
		1.5571281669726114 }, { make_pair(4, 5), 0.8732623126994838 }, {
		make_pair(4, 6), 1.108889969228049 }, { make_pair(4, 7),
		1.490954044606796 }, { make_pair(4, 8), 1.5505745831067692 }, {
		make_pair(4, 9), 1.6643412085537332 }, { make_pair(4, 10),
		1.521364594072136 }, { make_pair(4, 11), 1.3165336460366885 }, {
		make_pair(4, 12), 1.2989572026583633 }, { make_pair(4, 13),
		1.7376033486247136 }, { make_pair(4, 14), 1.0946271222686195 }, {
		make_pair(4, 15), 1.4093425168533116 }, { make_pair(4, 16),
		1.4412799026987895 }, { make_pair(4, 17), 1.7765038887875155 }, {
		make_pair(4, 18), 1.2745498917782174 }, { make_pair(4, 19),
		1.7629219871164792 }, { make_pair(5, 6), 1.289346236475462 }, {
		make_pair(5, 7), 1.1087882228915111 }, { make_pair(5, 8),
		1.5386383236095162 }, { make_pair(5, 9), 1.5864847448353807 }, {
		make_pair(5, 10), 1.4750214588048065 }, { make_pair(5, 11),
		1.2461612360696215 }, { make_pair(5, 12), 1.2607540372154786 }, {
		make_pair(5, 13), 1.6843769754160356 }, { make_pair(5, 14),
		1.1379547068550673 }, { make_pair(5, 15), 1.244344072759497 }, {
		make_pair(5, 16), 1.6352896598392934 }, { make_pair(5, 17),
		1.4381962213872233 }, { make_pair(5, 18), 1.5782891726771835 }, {
		make_pair(5, 19), 1.5126399522035463 }, { make_pair(6, 7),
		0.8733571175923197 }, { make_pair(6, 8), 1.587148486608941 }, {
		make_pair(6, 9), 1.5410748454710692 }, { make_pair(6, 10),
		1.2457316766898954 }, { make_pair(6, 11), 1.4767520675388714 }, {
		make_pair(6, 12), 1.6858462880724256 }, { make_pair(6, 13),
		1.2623340046978209 }, { make_pair(6, 14), 1.245141652888036 }, {
		make_pair(6, 15), 1.1390777466756155 }, { make_pair(6, 16),
		1.5577440060897845 }, { make_pair(6, 17), 1.697355511680117 }, {
		make_pair(6, 18), 1.2960156060995955 }, { make_pair(6, 19),
		1.3521004943160027 }, { make_pair(7, 8), 1.6630160105250575 }, {
		make_pair(7, 9), 1.551139320666103 }, { make_pair(7, 10),
		1.3162693574744762 }, { make_pair(7, 11), 1.5231191084036602 }, {
		make_pair(7, 12), 1.736953277259215 }, { make_pair(7, 13),
		1.2982404980919553 }, { make_pair(7, 14), 1.4094464587123365 }, {
		make_pair(7, 15), 1.0948515349150703 }, { make_pair(7, 16),
		1.8083982566714845 }, { make_pair(7, 17), 1.4207097378439073 }, {
		make_pair(7, 18), 1.6504164338651262 }, { make_pair(7, 19),
		1.105980844578514 }, { make_pair(8, 9), 1.2374365479449447 }, {
		make_pair(8, 10), 1.1987711448366474 }, { make_pair(8, 11),
		1.6941928740227394 }, { make_pair(8, 12), 1.24792187933318 }, {
		make_pair(8, 13), 0.7544534351220484 }, { make_pair(8, 14),
		1.2084059560581129 }, { make_pair(8, 15), 1.6252550146089235 }, {
		make_pair(8, 16), 1.5101185389321503 }, { make_pair(8, 17),
		1.622604351345613 }, { make_pair(8, 18), 1.495319331155472 }, {
		make_pair(8, 19), 1.4681979131184522 }, { make_pair(9, 10),
		1.6948412872717478 }, { make_pair(9, 11), 1.1985649615001661 }, {
		make_pair(9, 12), 0.7544150299805361 }, { make_pair(9, 13),
		1.248248873213137 }, { make_pair(9, 14), 1.6251189077959203 }, {
		make_pair(9, 15), 1.2095265891143607 }, { make_pair(9, 16),
		1.5890615846292173 }, { make_pair(9, 17), 1.6652858162445938 }, {
		make_pair(9, 18), 1.360959048297134 }, { make_pair(9, 19),
		1.3112745960332362 }, { make_pair(10, 11), 1.316002523695425 }, {
		make_pair(10, 12), 1.6974600303763911 }, { make_pair(10, 13),
		1.2095796762446704 }, { make_pair(10, 14), 0.7875916814048892 }, {
		make_pair(10, 15), 1.214038068389118 }, { make_pair(10, 16),
		1.5029561163419591 }, { make_pair(10, 17), 1.5909913198757268 }, {
		make_pair(10, 18), 1.6986349770052105 }, { make_pair(10, 19),
		1.4961322180216545 }, { make_pair(11, 12), 1.2085056067493043 }, {
		make_pair(11, 13), 1.6951527898936165 }, { make_pair(11, 14),
		1.2143622250345927 }, { make_pair(11, 15), 0.7879118215850902 }, {
		make_pair(11, 16), 1.5263933439051656 }, { make_pair(11, 17),
		1.6352297920980154 }, { make_pair(11, 18), 1.5948775504089308 }, {
		make_pair(11, 19), 1.5319860984358884 }, { make_pair(12, 13),
		1.437068750361695 }, { make_pair(12, 14), 1.694794469239858 }, {
		make_pair(12, 15), 1.1702728833431355 }, { make_pair(12, 16),
		1.5573897315418705 }, { make_pair(12, 17), 1.7850444354851502 }, {
		make_pair(12, 18), 1.4359919035052873 }, { make_pair(12, 19),
		1.595563410663024 }, { make_pair(13, 14), 1.1700367366804003 }, {
		make_pair(13, 15), 1.6957534476741514 }, { make_pair(13, 16),
		1.6774090100030388 }, { make_pair(13, 17), 1.583794018728718 }, {
		make_pair(13, 18), 1.5823747262305339 }, { make_pair(13, 19),
		1.3085850503357264 }, { make_pair(14, 15), 1.5340642704445115 }, {
		make_pair(14, 16), 1.4726438094417185 }, { make_pair(14, 17),
		1.4782852692640378 }, { make_pair(14, 18), 1.661858231287885 }, {
		make_pair(14, 19), 1.5810982937458442 }, { make_pair(15, 16),
		1.5478113442579584 }, { make_pair(15, 17), 1.6758864884380875 }, {
		make_pair(15, 18), 1.5057337091627314 }, { make_pair(15, 19),
		1.3399915570941325 }, { make_pair(16, 17), 1.217427036513581 }, {
		make_pair(16, 18), 1.331202341456014 }, { make_pair(16, 19),
		1.4655979135540935 }, { make_pair(17, 18), 1.83400283089752 }, {
		make_pair(17, 19), 1.319848059810783 }, { make_pair(18, 19),
		1.2147744180018427 } };

// need to add more tests once we can assign types properly in SDF files.

/**
 * Test rms values for a multiconformer file, using the specialized multiconformer
 * evaluator.  Reference values from java program.
 */
TEST_F(LeastSquaresTest, TestMultiConformerOmegaFitting) {

	const string inFile = "../resources/1g9v_omega_20.mol2";
	auto mols = MulticonformerMolecule::readMoleculesFromFile(inFile);
	ASSERT_EQ(mols.size(), static_cast<size_t>(1));
	auto nConformers = mols[0]->nConformers();
	ASSERT_EQ(nConformers, static_cast<size_t>(20));
	MultiConformerRmsEvalution rmsEvaluation(*mols[0]);
	const auto & rmsArray = rmsEvaluation.evaluateRms();

	for (size_t i = 0; i < nConformers; i++) {
		for (size_t j = i + 1; j < nConformers; j++) {
			const auto rmsd = rmsArray(i, j);
			auto ref = omegaResults.at(make_pair(i, j));
			REPORT(Reporter::NORMAL) << "Rmsd between conformers " << i + 1
					<< " and " << j + 1 << " " << rmsd;
			ASSERT_TRUE(equals(ref, rmsd, TOLERANCE));
		}
	}

	// build nearest neighbour lists and check they are good
	rmsEvaluation.determineNearestNeighbourLists();
	ASSERT_TRUE(rmsEvaluation.checkNeighbourLists());
}

/**
 * Test rms values for multiconformer files.  Reference values from java program.
 */
TEST_F(LeastSquaresTest, TestLigandToMultiConformer) {

	// reference values obtained using: java com.arenapharm.molecule.FindClosestConformer ../resources/1gkc_ligand_charged.mol ../resources/1gkc_omega.mol2

	// First treat multiconformer file as series of molecules
	auto ligandMolecule = Molecule::readAndInitializeMolecule(
			"../resources/1gkc_ligand_charged.mol");
	auto conformers = Molecule::readMoleculesFromFile(
			"../resources/1gkc_omega.mol2");
	auto minRms = numeric_limits<double>::max();
	for (const auto & conformer : conformers) {
		MolecularRms molecularRms(*ligandMolecule, *conformer);

		molecularRms.determineRms();
		auto rmsd = molecularRms.getMinRmsd();
		REPORT(Reporter::NORMAL) << "Rmsd between ligand and conformer "
				<< rmsd;
		if (rmsd < minRms) {
			minRms = rmsd;
		}
	}

	REPORT(Reporter::NORMAL) << "Minimim rmsd between ligand and conformer "
			<< minRms;
	ASSERT_TRUE(equals(minRms, 0.38580454355102456, 5));

	// the same thing reading multiconformer molecules
	minRms = numeric_limits<double>::max();
	auto multiMolecule = MulticonformerMolecule::readAndInitializeMolecule(
			"../resources/1gkc_omega.mol2");
	MolecularRms molecularRms(*ligandMolecule, *multiMolecule);
	// only determine isomorphisms once
	molecularRms.determineIsomorphisms();
	for (const auto & conformer : multiMolecule->getConformers()) {
		auto rmsd = molecularRms.determineRms(ligandMolecule->getCoords(),
				conformer->getCoordinates());
		REPORT(Reporter::NORMAL) << "Rmsd between ligand and conformer "
				<< rmsd;
		if (rmsd < minRms) {
			minRms = rmsd;
		}
	}

	REPORT(Reporter::NORMAL) << "Minimim rmsd between ligand and conformer "
			<< minRms;
	ASSERT_TRUE(equals(minRms, 0.38580454355102456, 5));
}
/**
 * Test rms values for multiconformer files.  Reference values from java program.
 */
TEST_F(LeastSquaresTest, TestOmegaFitting) {

	// reference values obtained using: java com.arenapharm.molecule.EvaluateConformerRms ~/src/gape/resources/1g9v_omega_20.sdf

	const string inFile = "../resources/1g9v_omega_20.mol2";
	auto molecules = Molecule::readMoleculesFromFile(inFile);
	auto size = molecules.size();
	ASSERT_EQ(size, static_cast<size_t>(20));
	for (size_t i = 0; i < size; i++) {
		for (size_t j = i + 1; j < size; j++) {
			MolecularRms molecularRms(*molecules.at(j), *molecules.at(i));
			molecularRms.determineRms();
			auto rmsd = molecularRms.getMinRmsd();
			REPORT(Reporter::NORMAL) << "Rmsd between conformers " << i + 1
					<< " and " << j + 1 << " " << rmsd;
			REPORT(Reporter::NORMAL) << "All rmsds "
					<< collectionToString(molecularRms.getIsomorphismRmsds());
			auto ref = omegaResults.at(make_pair(i, j));
			ASSERT_TRUE(equals(ref, rmsd, TOLERANCE));
		}
	}

}

/**
 * Test of 1bpd1bqm self-fitting- useful for debugging
 */
TEST_F(LeastSquaresTest, TestBqmFitting) {

	const string inFile = "../resources/pdb1bqm.mol2";
	auto molecules = Molecule::readMoleculesFromFile(inFile);
	auto & mol = molecules.at(0);

	auto molecularRms = GarethMol::MolecularRms::createSelfComparer(*mol);
	molecularRms->determineRms();
	ASSERT_TRUE(equals(molecularRms->getMinRmsd(), .0, TOLERANCE));
	auto nZero = 0;
	for (auto rms : molecularRms->getIsomorphismRmsds()) {
		if (equals(rms, .0, TOLERANCE))
			nZero++;
	}
	ASSERT_EQ(nZero, 1);

}

/**
 * Test of some hiv inhibitors.
 */
TEST_F(LeastSquaresTest, CompareHiv) {
	const string inFile = "../resources/hiv_rt2.mol2";
	auto molecules = Molecule::readMoleculesFromFile(inFile);
	for (auto & mol : molecules) {
		REPORT(Reporter::NORMAL) << "Processing molecule " << mol->getName();
		auto molecularRms = GarethMol::MolecularRms::createSelfComparer(*mol);
		molecularRms->determineRms();
		auto rmsd = molecularRms->getMinRmsd();
		REPORT(Reporter::NORMAL) << "Min rmsd for " << mol->getName() << " is "
				<< rmsd;
		ASSERT_TRUE(equals(rmsd, .0, TOLERANCE));
		auto nZero = 0;
		for (auto rms : molecularRms->getIsomorphismRmsds()) {
			if (equals(rms, .0, TOLERANCE))
				nZero++;
		}
		ASSERT_EQ(nZero, 1);
	}
}

/**
 * Test of benzene self-fitting
 */
TEST_F(LeastSquaresTest, TestBezeneFitting) {

	const string inFile = "../resources/benzene.mol2";
	auto molecules = Molecule::readMoleculesFromFile(inFile);
	auto & mol = molecules.at(0);

	auto molecularRms = GarethMol::MolecularRms::createSelfComparer(*mol);
	// for benzene all 12 isomorphisms should have an rms of 0
	molecularRms->determineRms();
	ASSERT_TRUE(equals(molecularRms->getMinRmsd(), .0, TOLERANCE));
	for (auto rms : molecularRms->getIsomorphismRmsds()) {
		ASSERT_TRUE(equals(rms, .0, TOLERANCE));
	}
}

}
