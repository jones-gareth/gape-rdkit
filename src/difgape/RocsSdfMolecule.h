/*
 * RocsSdfMolecule.h
 *
 *  Created on: May 14, 2014
 *      Author: Gareth Jones
 *
 *  Used to store a summary of an Openeye molecule from a ROCS run.
 *  This class identifies features for use in defining distance constraints.
 *
 *  I've been unable to get the serialization stuff working yet
 */

#ifndef ROCSSDFMOLECULE_H_
#define ROCSSDFMOLECULE_H_

/*
 couldn't get boost serialization to work.
 #include <boost/serialization/access.hpp>
 #include <boost/serialization/string.hpp>
 #include <boost/serialization/array.hpp>
 #include <boost/serialization/vector.hpp>
 #include <boost/serialization/utility.hpp>
 #include <boost/serialization/map.hpp>
 #include <boost/serialization/base_object.hpp>
 #include <boost/serialization/list.hpp>
 #include <boost/serialization/shared_ptr.hpp>
 #include <boost/archive/text_iarchive.hpp>
 #include <boost/archive/text_oarchive.hpp>
 #include <boost/smart_ptr/shared_ptr.hpp>
 */
#include <cassert>
#include <map>
#include <memory>
#include <string>
#include <vector>
#include <array>

#include "../mol/mol.h"
#include "../mol/Molecule.h"
// Don't include RocsSdfParser parser here we get into nested dependencies issues
//#include "RocsSdfParser.h"

namespace Difgape {

using namespace std;
using namespace GarethMol;

enum class ScoreName
;

class RocsSdfMolecule {
public:
	//using RocsSdfMoleculePtr = boost::shared_ptr<RocsSdfMolecule>;
	using RocsSdfMoleculePtr = shared_ptr<RocsSdfMolecule>;

	/**
	 * This constructor to be used for serialization only
	 */
	RocsSdfMolecule() {
	}

	/**
	 * Create the summary molecule from a normal molecule.
	 * @param molecule
	 */
	RocsSdfMolecule(const Molecule::MoleculePtr & molecule);

	virtual ~RocsSdfMolecule() {
	}
	;

	/**
	 * @return the conformer number
	 */
	const int getConformerNo() const {
		return conformerNo;
	}

	const size_t getNFeatures() const {
		return featureCoordinates.size();
	}
	/**
	 *
	 * @return feature coordinates
	 */
	const CoordMatrix & getFeatureCoordinates() const {
		return featureCoordinates;
	}

	/**
	 *
	 * @return feature names
	 */
	const vector<string>& getFeatureNames() const {
		return featureNames;
	}

	const map<ScoreName, double>& getScores() const {
		return scores;
	}

	/**
	 *
	 * @return structure name
	 */
	const string& getStructureName() const {
		return structureName;
	}

	const string getName() const {
		string name = structureName + "_" + to_string(conformerNo);
		return name;
	}
	/**
	 *
	 * @param no
	 * @return coordinates for  a particular feature
	 */
	const CoordVector getFeatureCoordinate(int no) const {
		return featureCoordinates.col(no);
	}

	/**
	 *
	 * @param no
	 * @return the name for a particular feature
	 */
	const string & getFeatureName(int no) const {
		return featureNames.at(no);
	}

	/**
	 * Sets the conformer no. This is required as the new ROCS/Omega does not
	 * include conformer no in initial query.
	 *
	 * @param conformerNo
	 */
	void setConformerNo(int conformerNo) {
		this->conformerNo = conformerNo;
	}

	/**
	 * sets a particular score
	 *
	 * @param scoreName
	 * @param score
	 */
	void setScore(Difgape::ScoreName scoreName, double score);

	/**
	 *
	 * @param scoreName
	 * @return a particular score
	 */
	const double getScore(ScoreName scoreName) const {
		return scores.at(scoreName);
	}
	/**
	 *
	 * @return the number of features
	 */
	const int noFeatures() const {
		return featureNames.size();
	}

	/**
	 * Check that this conformer has equivalent features to the reference
	 * strucure.
	 *
	 * @param other
	 * @return
	 */
	bool checkReference(RocsSdfMoleculePtr & other) const;

private:

	// serialization- wont work
	/*
	 friend class boost::serialization::access;
	 template<class Archive>
	 void serialize(Archive & ar, const unsigned int version) {
	 ar & conformerNo;
	 ar & structureName;
	 ar & featureCoordinates;
	 ar & featureNames;
	 ar & scores;
	 }
	 */

	string structureName { "" };
	int conformerNo { 0 };
	CoordMatrix featureCoordinates { };
	vector<string> featureNames { };
	map<ScoreName, double> scores;

	/**
	 * Find features and feature coordinates for distance constraints. Currently
	 * use hetero atoms. It would be best to identify rigid portions of
	 * molecules and use centroids.
	 *
	 * @param molecule
	 */
	void findFeatures(const Molecule::MoleculePtr & molecule);

	RocsSdfMolecule(const RocsSdfMolecule &) = delete;
	RocsSdfMolecule & operator =(const RocsSdfMolecule &) = delete;
};

} /* namespace Difgape */

// allow std::array to be serialized
/*
 namespace boost {
 namespace serialization {

 template<class Archive, class T, size_t N>
 void serialize(Archive & ar, std::array<T, N> & a, const unsigned int version) {
 ar & boost::serialization::make_array(a.data(), a.size());
 }

 } // namespace serialization
 } // namespace boost
 */
#endif /* ROCSSDFMOLECULE_H_ */
