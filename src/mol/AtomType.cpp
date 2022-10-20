/*
 * AtomType.cpp
 *
 *  Created on: Jun 21, 2013
 *      Author: Gareth Jones
 */

#include "AtomType.h"
#include "Reporter.h"
#include <algorithm>
#include <utility>
#include <boost/format.hpp>

namespace GarethMol {

using namespace GarethUtil;
using namespace std;

AtomType::~AtomType() {
}

/**
 * Determine the gaussian width for shape overlay from a radius
 *
 * @param radius
 * @return
 */
double AtomType::determineGaussianWidth(const double radius) {
    double r = radius;
    double alpha = pow((3.0 * ATOMIC_N) / (4 * M_PI * r * r * r),
            2.0 / 3.0) * M_PI;
    return alpha;
}

struct AtomType::StaticData AtomType::buildTypes() {

    TypesVector typesVector;
    typesVector.reserve(100);
    //TODO should define all these guys in an external file- see charge class
    typesVector.push_back(
            make_unique<const AtomType>(AtomTypeId::C3, "C.3", 1.7,
                    Geometry::TET, 12.0, vector<int> { 4 }, 6));
    typesVector.push_back(
            make_unique<const AtomType>(AtomTypeId::C2, "C.2", 1.7,
                    Geometry::TRI, 12.0, vector<int> { 4 }, 6));
    typesVector.push_back(
            make_unique<const AtomType>(AtomTypeId::C1, "C.1", 1.7,
                    Geometry::LIN, 12.0, vector<int> { 4 }, 6));
    typesVector.push_back(
            make_unique<const AtomType>(AtomTypeId::CAR, "C.ar", 1.7,
                    Geometry::TRI, 12.0, vector<int> { 4 }, 6));
    typesVector.push_back(
            make_unique<const AtomType>(AtomTypeId::CCAT, "C.cat", 1.7,
                    Geometry::TRI, 12.0, vector<int> { 4 }, 6));
    typesVector.push_back(
            make_unique<const AtomType>(AtomTypeId::C, "C", 1.7, Geometry::NONE,
                    12.0, vector<int> { 4 }, 6));
    typesVector.push_back(
            make_unique<const AtomType>(AtomTypeId::N3, "N.3", 1.55,
                    Geometry::TET, 14.0, vector<int> { 3 }, 7));
    typesVector.push_back(
            make_unique<const AtomType>(AtomTypeId::N2, "N.2", 1.55,
                    Geometry::TRI, 14.0, vector<int> { 3 }, 7));
    typesVector.push_back(
            make_unique<const AtomType>(AtomTypeId::N1, "N.1", 1.55,
                    Geometry::LIN, 14.0, vector<int> { 3 }, 7));
    typesVector.push_back(
            make_unique<const AtomType>(AtomTypeId::NAR, "N.ar", 1.55,
                    Geometry::TRI, 14.0, vector<int> { 3 }, 7));
    typesVector.push_back(
            make_unique<const AtomType>(AtomTypeId::NAM, "N.am", 1.55,
                    Geometry::TRI, 14.0, vector<int> { 3 }, 7));
    typesVector.push_back(
            make_unique<const AtomType>(AtomTypeId::NPL3, "N.pl3", 1.55,
                    Geometry::TRI, 14.0, vector<int> { 3 }, 7));
    typesVector.push_back(
            make_unique<const AtomType>(AtomTypeId::N4, "N.4", 1.55,
                    Geometry::TET, 14.0, vector<int> { 3 }, 7));
    typesVector.push_back(
            make_unique<const AtomType>(AtomTypeId::N, "N", 1.55,
                    Geometry::NONE, 14.0, vector<int> { 3 }, 7));
    typesVector.push_back(
            make_unique<const AtomType>(AtomTypeId::O3, "O.3", 1.52,
                    Geometry::TET, 16.0, vector<int> { 2 }, 8));
    typesVector.push_back(
            make_unique<const AtomType>(AtomTypeId::O2, "O.2", 1.52,
                    Geometry::TRI, 16.0, vector<int> { 2 }, 8));
    typesVector.push_back(
            make_unique<const AtomType>(AtomTypeId::OCO2, "O.co2", 1.52,
                    Geometry::TRI, 16.0, vector<int> { 2 }, 8));
    typesVector.push_back(
            make_unique<const AtomType>(AtomTypeId::OAR, "O.ar", 1.52,
                    Geometry::TRI, 16.0, vector<int> { 2 }, 8));
    typesVector.push_back(
            make_unique<const AtomType>(AtomTypeId::OANY, "O", 1.52,
                    Geometry::NONE, 16.0, vector<int> { 2 }, 8));
    typesVector.push_back(
            make_unique<const AtomType>(AtomTypeId::S3, "S.3", 1.8,
                    Geometry::TET, 32.0, vector<int> { 2, 4, 6 }, 16));
    typesVector.push_back(
            make_unique<const AtomType>(AtomTypeId::S2, "S.2", 1.8,
                    Geometry::TRI, 32.0, vector<int> { 2, 4, 6 }, 16));
    typesVector.push_back(
            make_unique<const AtomType>(AtomTypeId::SO, "S.o", 1.7,
                    Geometry::TET, 32.0, vector<int> { 2, 4, 6 }, 16));
    typesVector.push_back(
            make_unique<const AtomType>(AtomTypeId::SO2, "S.o2", 1.7,
                    Geometry::TET, 32.0, vector<int> { 2, 4, 6 }, 16));
    typesVector.push_back(
            make_unique<const AtomType>(AtomTypeId::S, "S", 1.7, Geometry::NONE,
                    32.0, vector<int> { 2, 4, 6 }, 16));
    typesVector.push_back(
            make_unique<const AtomType>(AtomTypeId::P3, "P.3", 1.8,
                    Geometry::TET, 32.0, vector<int> { 3, 5 }, 15));
    typesVector.push_back(
            make_unique<const AtomType>(AtomTypeId::P, "P", 1.8, Geometry::NONE,
                    32.0, vector<int> { 3, 5 }, 15));
    typesVector.push_back(
            make_unique<const AtomType>(AtomTypeId::H, "H", 1.5, Geometry::LIN,
                    1.0, vector<int> { 1 }, 1));
    typesVector.push_back(
            make_unique<const AtomType>(AtomTypeId::CL, "Cl", 1.75,
                    Geometry::TET, 34.4, vector<int> { 1 }, 17));
    typesVector.push_back(
            make_unique<const AtomType>(AtomTypeId::BR, "Br", 1.85,
                    Geometry::TET, 79.9, vector<int> { 1 }, 35));
    typesVector.push_back(
            make_unique<const AtomType>(AtomTypeId::I, "I", 1.98, Geometry::TET,
                    126.9, vector<int> { 1 }, 53));
    typesVector.push_back(
            make_unique<const AtomType>(AtomTypeId::SI, "Si", 1.2,
                    Geometry::TET, 28.1, vector<int> { 4 }, 14));
    typesVector.push_back(
            make_unique<const AtomType>(AtomTypeId::LP, "Lp", 1.0,
                    Geometry::LIN, 0, vector<int> { 1 }, 0));
    typesVector.push_back(
            make_unique<const AtomType>(AtomTypeId::DU, "Du", 1.0,
                    Geometry::NONE, 0, vector<int> { }, 0));
    typesVector.push_back(
            make_unique<const AtomType>(AtomTypeId::NA, "Na", 1.2, Geometry::CO,
                    23.0, vector<int> { }, 11));
    typesVector.push_back(
            make_unique<const AtomType>(AtomTypeId::K, "K", 1.2, Geometry::CO,
                    39.1, vector<int> { }, 19));
    typesVector.push_back(
            make_unique<const AtomType>(AtomTypeId::F, "F", 1.47, Geometry::TET,
                    19.0, vector<int> { 1 }, 9));
    typesVector.push_back(
            make_unique<const AtomType>(AtomTypeId::CA, "Ca", 1.11,
                    Geometry::CO, 40.1, vector<int> { }, 20));
    typesVector.push_back(
            make_unique<const AtomType>(AtomTypeId::LI, "Li", 1.2, Geometry::CO,
                    6.94, vector<int> { }, 3));
    typesVector.push_back(
            make_unique<const AtomType>(AtomTypeId::AL, "Al", 1.2, Geometry::CO,
                    27.0, vector<int> { }, 13));
    typesVector.push_back(
            make_unique<const AtomType>(AtomTypeId::MG, "Mg", 0.74,
                    Geometry::CO, 24.3, vector<int> { }, 12));
    typesVector.push_back(
            make_unique<const AtomType>(AtomTypeId::ZN, "Zn", 0.83,
                    Geometry::CO, 65.4, vector<int> { }, 30));
    typesVector.push_back(
            make_unique<const AtomType>(AtomTypeId::FE, "Fe", 0.83,
                    Geometry::CO, 55.8, vector<int> { }, 26));
    typesVector.push_back(
            make_unique<const AtomType>(AtomTypeId::MN, "Mn", 0.89,
                    Geometry::CO, 54.9, vector<int> { }, 25));
    typesVector.push_back(
            make_unique<const AtomType>(AtomTypeId::B, "B", 1.70, Geometry::CO,
                    10.81, vector<int> { }, 5));
    typesVector.push_back(
            make_unique<const AtomType>(AtomTypeId::HET, "Het", 0.0,
                    Geometry::TET, 0.0, vector<int> { }, 0));
    typesVector.push_back(
            make_unique<const AtomType>(AtomTypeId::HAL, "Hal", 0.0,
                    Geometry::NONE, 0.0, vector<int> { }, 0));
    typesVector.push_back(
            make_unique<const AtomType>(AtomTypeId::HEV, "Hev", 0.0,
                    Geometry::NONE, 0.0, vector<int> { }, 0));
    typesVector.push_back(
            make_unique<const AtomType>(AtomTypeId::DUC, "Du.C", 0.0,
                    Geometry::NONE, 0.0, vector<int> { }, 0));
    typesVector.push_back(
            make_unique<const AtomType>(AtomTypeId::ANY, "Any", 0.0,
                    Geometry::NONE, 0.0, vector<int> { }, 0));
    typesVector.push_back(
            make_unique<const AtomType>(AtomTypeId::WILD, "*", 0.0,
                    Geometry::NONE, 0.0, vector<int> { }, 0));
    typesVector.push_back(
            make_unique<const AtomType>(AtomTypeId::ATM_NONE, "None", 0.0,
                    Geometry::NONE, 0.0, vector<int> { }, 0));
    typesVector.push_back(
            make_unique<const AtomType>(AtomTypeId::CU, "Cu", 1.40,
                    Geometry::CO, 65.3, vector<int> { }, 29));
    typesVector.push_back(
            make_unique<const AtomType>(AtomTypeId::SE, "Se", 1.90,
                    Geometry::CO, 78.9, vector<int> { }, 34));
    typesVector.push_back(
            make_unique<const AtomType>(AtomTypeId::PB, "Pb", 2.02,
                    Geometry::CO, 207.2, vector<int> { }, 82));
    typesVector.push_back(
            make_unique<const AtomType>(AtomTypeId::HG, "Hg", 1.55,
                    Geometry::CO, 200.6, vector<int> { }, 80));
    typesVector.push_back(
            make_unique<const AtomType>(AtomTypeId::SB, "Sb", -1, Geometry::CO,
                    121.8, vector<int> { }, 51));
    typesVector.push_back(
            make_unique<const AtomType>(AtomTypeId::BI, "Bi", -1, Geometry::CO,
                    209.0, vector<int> { }, 83));
    typesVector.push_back(
            make_unique<const AtomType>(AtomTypeId::TI, "Ti", -1, Geometry::CO,
                    47.8, vector<int> { }, 22));
    typesVector.push_back(
            make_unique<const AtomType>(AtomTypeId::AS, "As", 1.85,
                    Geometry::CO, 74.9, vector<int> { }, 22));
    typesVector.push_back(
            make_unique<const AtomType>(AtomTypeId::NI, "Ni", 1.63,
                    Geometry::CO, 58.7, vector<int> { }, 28));
    typesVector.push_back(
            make_unique<const AtomType>(AtomTypeId::BA, "Ba", 1.63,
                    Geometry::CO, 137.3, vector<int> { }, 56));
    typesVector.push_back(
            make_unique<const AtomType>(AtomTypeId::AG, "Ag", 1.63,
                    Geometry::CO, 107.9, vector<int> { }, 47));
    typesVector.push_back(
            make_unique<const AtomType>(AtomTypeId::SN, "Sn", 1.90,
                    Geometry::CO, 118.7, vector<int> { }, 50));
    typesVector.push_back(
            make_unique<const AtomType>(AtomTypeId::CR, "Cr", -1, Geometry::CO,
                    52.0, vector<int> { }, 24));
    typesVector.push_back(
            make_unique<const AtomType>(AtomTypeId::CD, "Cd", 1.58,
                    Geometry::CO, 112.4, vector<int> { }, 48));
    typesVector.push_back(
            make_unique<const AtomType>(AtomTypeId::CO, "Co", -1, Geometry::CO,
                    59.0, vector<int> { }, 27));
    typesVector.push_back(
            make_unique<const AtomType>(AtomTypeId::AU, "Au", 1.66,
                    Geometry::CO, 197.0, vector<int> { }, 79));
    typesVector.push_back(
            make_unique<const AtomType>(AtomTypeId::CE, "Ce", -1, Geometry::CO,
                    140.1, vector<int> { }, 58));
    typesVector.push_back(
            make_unique<const AtomType>(AtomTypeId::BE, "Be", -1, Geometry::CO,
                    9.01, vector<int> { }, 4));
    typesVector.push_back(
            make_unique<const AtomType>(AtomTypeId::TH, "Th", -1, Geometry::CO,
                    232.0, vector<int> { }, 90));
    typesVector.push_back(
            make_unique<const AtomType>(AtomTypeId::V, "V", -1, Geometry::CO,
                    50.9, vector<int> { }, 23));
    typesVector.push_back(
            make_unique<const AtomType>(AtomTypeId::ZR, "Zr", -1, Geometry::CO,
                    91.2, vector<int> { }, 40));
    typesVector.push_back(
            make_unique<const AtomType>(AtomTypeId::PT, "Pt", 1.75,
                    Geometry::CO, 195.1, vector<int> { }, 78));
    typesVector.push_back(
            make_unique<const AtomType>(AtomTypeId::TE, "Te", 2.06,
                    Geometry::CO, 127.6, vector<int> { }, 52));
    typesVector.push_back(
            make_unique<const AtomType>(AtomTypeId::RA, "Ra", -1, Geometry::CO,
                    226.0, vector<int> { }, 88));
    typesVector.push_back(
            make_unique<const AtomType>(AtomTypeId::GE, "Ge", -1, Geometry::CO,
                    72.6, vector<int> { }, 32));
    typesVector.push_back(
            make_unique<const AtomType>(AtomTypeId::GA, "Ga", 1.87,
                    Geometry::CO, 69.7, vector<int> { }, 31));
    typesVector.push_back(
            make_unique<const AtomType>(AtomTypeId::TL, "Tl", 1.96,
                    Geometry::CO, 204.4, vector<int> { }, 81));
    typesVector.push_back(
            make_unique<const AtomType>(AtomTypeId::CS, "Cs", -1, Geometry::CO,
                    132.9, vector<int> { }, 55));
    typesVector.push_back(
            make_unique<const AtomType>(AtomTypeId::IN, "In", 1.93,
                    Geometry::CO, 114.8, vector<int> { }, 49));
    typesVector.push_back(
            make_unique<const AtomType>(AtomTypeId::SM, "Sm", -1, Geometry::CO,
                    150.3, vector<int> { }, 62));
    typesVector.push_back(
            make_unique<const AtomType>(AtomTypeId::PD, "Pd", 1.63,
                    Geometry::CO, 106.4, vector<int> { }, 46));

    std::map<const AtomType::AtomTypeId, const AtomType *> typesLookup;
    std::map<const std::string, const AtomType *> typesFromNameLookup;

    for (const auto & atomType : typesVector) {
        typesLookup.insert(
                pair<AtomType::AtomTypeId, const AtomType *>(
                        atomType->getType(), atomType.get()));
        typesFromNameLookup.insert(
                pair<string, const AtomType *>(atomType->getName(),
                        atomType.get()));
    }

    StaticData staticData { move(typesVector), typesLookup, typesFromNameLookup };
    return staticData;
}

AtomType::StaticData & AtomType::atomTypes() {
    static auto staticData = buildTypes();
    return staticData;
}

const std::vector<AtomType::AtomTypeId> AtomType::nonAtomTypes { AtomTypeId::X,
        AtomTypeId::Y, AtomTypeId::Z, AtomTypeId::DU, AtomTypeId::WILD,
        AtomTypeId::DUC, AtomTypeId::LP };

const AtomType & AtomType::typeFromName(const std::string & name) {

    // need these to read TAFF files
    if (equals(name, "HEV")) {
        return typeFromName("Hev");
    }
    if (equals(name, "LP")) {
        return typeFromName("Lp");
    }
    if (equals(name, "ANY")) {
        return typeFromName("Any");
    }
    if (equals(name, "HAL")) {
        return typeFromName("Hal");
    }
    if (equals(name, "HET")) {
        return typeFromName("Het");
    }

    const auto & typesFromNameLookup = atomTypes().typesFromNameLookup;
    auto iter = typesFromNameLookup.find(name);
    if (iter == typesFromNameLookup.end()) {
        REPORT(Reporter::WARN) << "No type for Sybyl atom " << name;
        return typeFromName("Du");
    } else {
        return *iter->second;
    }
}

const std::string & AtomType::nameFromTypeId(const AtomTypeId & type) {
    const auto & typesLookup = atomTypes().typesLookup;
    assert(typesLookup.count(type) == 1);
    return typesLookup.at(type)->getName();

}

const AtomType & AtomType::typeFromTypeId(const AtomTypeId & type) {
    const auto & typesLookup = atomTypes().typesLookup;
    return *typesLookup.at(type);
}

const bool AtomType::isOxygenType(const AtomTypeId & type) {
    if (type == AtomTypeId::O2 || type == AtomTypeId::O3
            || type == AtomTypeId::OCO2 || type == AtomTypeId::OANY
            || type == AtomTypeId::OAR)
        return true;
    return false;
}

const bool AtomType::isHeteroType(const AtomTypeId & type) {
    if (isOxygenType(type) || isNitrogenType(type) || isPhosphorousType(type)
            || isSulphurType(type))
        return true;
    if (isHalogen(type))
        return true;
    return false;
}

const bool AtomType::isHalogen(const AtomTypeId & type) {
    if (type == AtomTypeId::CL || type == AtomTypeId::BR
            || type == AtomTypeId::F)
        return true;
    return false;
}

const bool AtomType::isCarbonType(const AtomTypeId & type) {
    if (type == AtomTypeId::C1 || type == AtomTypeId::C2
            || type == AtomTypeId::C3 || type == AtomTypeId::CAR
            || type == AtomTypeId::CCAT || type == AtomTypeId::C)
        return true;
    return false;
}

const bool AtomType::isSulphurType(const AtomTypeId & type) {
    if (type == AtomTypeId::S2 || type == AtomTypeId::S3
            || type == AtomTypeId::SO || type == AtomTypeId::SO2
            || type == AtomTypeId::S)
        return true;
    return false;
}

const bool AtomType::isPhosphorousType(const AtomTypeId & type) {
    if (type == AtomTypeId::P3 || type == AtomTypeId::P)
        return true;
    return false;
}

const bool AtomType::isNitrogenType(const AtomTypeId & type) {
    if (type == AtomTypeId::N1 || type == AtomTypeId::N2
            || type == AtomTypeId::N3 || type == AtomTypeId::NAR
            || type == AtomTypeId::N4 || type == AtomTypeId::NPL3
            || type == AtomTypeId::NAM || type == AtomTypeId::N)
        return true;
    return false;
}

const bool AtomType::isHeavy(const AtomTypeId & type) {
    return isRealAtom(type) && type != AtomTypeId::H;
}

const bool AtomType::isRealAtom(const AtomTypeId & type) {
    auto iter = std::find(nonAtomTypes.cbegin(), nonAtomTypes.cend(), type);
    return iter == nonAtomTypes.cend();
}

const bool AtomType::isNotDummy(const AtomTypeId & type) {
    if (type == AtomTypeId::DU || type == AtomTypeId::NA
            || type == AtomTypeId::K || type == AtomTypeId::CA
            || type == AtomTypeId::LI || type == AtomTypeId::MG
            || type == AtomTypeId::ZN || type == AtomTypeId::FE
            || type == AtomTypeId::MN || type == AtomTypeId::LP)
        return false;
    return true;
}

const bool AtomType::matchType(const AtomTypeId & type1,
        const AtomTypeId & type2, const bool ignoreAromatic) {

    if (type1 == AtomTypeId::ANY || type2 == AtomTypeId::ANY)
        return true;
    if (type1 == AtomTypeId::WILD || type2 == AtomTypeId::WILD)
        return true;

    if (type1 == AtomTypeId::X && isHeavy(type2))
        return true;
    if (type1 == AtomTypeId::Y && isHeavy(type2))
        return true;
    if (type1 == AtomTypeId::Z && isHeavy(type2))
        return true;
    if (type2 == AtomTypeId::X && isHeavy(type1))
        return true;
    if (type2 == AtomTypeId::Y && isHeavy(type1))
        return true;
    if (type2 == AtomTypeId::Z && isHeavy(type1))
        return true;

    if (type1 == AtomTypeId::OANY && isOxygenType(type2))
        return true;
    if (type2 == AtomTypeId::OANY && isOxygenType(type1))
        return true;
    if (type1 == AtomTypeId::C && isCarbonType(type2))
        return true;
    if (type2 == AtomTypeId::C && isCarbonType(type1))
        return true;
    if (type1 == AtomTypeId::N && isNitrogenType(type2))
        return true;
    if (type2 == AtomTypeId::N && isNitrogenType(type1))
        return true;
    if (type1 == AtomTypeId::P && isPhosphorousType(type2))
        return true;
    if (type2 == AtomTypeId::P && isPhosphorousType(type1))
        return true;
    if (type1 == AtomTypeId::S && isSulphurType(type2))
        return true;
    if (type2 == AtomTypeId::S && isSulphurType(type1))
        return true;

    if (type1 == AtomTypeId::HEV && isHeavy(type1))
        return true;
    if (type2 == AtomTypeId::HEV && isHeavy(type2))
        return true;

    if (type1 == AtomTypeId::HAL
            && (type2 == AtomTypeId::CL || type2 == AtomTypeId::BR
                    || type2 == AtomTypeId::F))
        return true;
    if (type2 == AtomTypeId::HAL
            && (type1 == AtomTypeId::CL || type1 == AtomTypeId::BR
                    || type1 == AtomTypeId::F))
        return true;

    if (ignoreAromatic) {
        if (type1 == AtomTypeId::O3 && type2 == AtomTypeId::OAR)
            return true;
        if (type1 == AtomTypeId::CAR && type2 == AtomTypeId::C2)
            return true;
        if (type1 == AtomTypeId::NAR
                && (type2 == AtomTypeId::N2 || type2 == AtomTypeId::NPL3))
            return true;
        if (type2 == AtomTypeId::O3 && type1 == AtomTypeId::OAR)
            return true;
        if (type2 == AtomTypeId::CAR && type1 == AtomTypeId::C2)
            return true;
        if (type2 == AtomTypeId::NAR
                && (type1 == AtomTypeId::N2 || type2 == AtomTypeId::NPL3))
            return true;
    }
    if (type1 == type2)
        return true;

    return false;
}

const bool AtomType::matchElementalType(const AtomTypeId & type1,
        const AtomTypeId & type2, bool ignoreAromatic) {
    if (isOxygenType(type1) && isOxygenType(type2))
        return true;
    if (isCarbonType(type1) && isCarbonType(type2))
        return true;
    if (isNitrogenType(type1) && isNitrogenType(type2))
        return true;
    if (isPhosphorousType(type2) && isPhosphorousType(type1))
        return true;
    if (isSulphurType(type1) && isSulphurType(type2))
        return true;

    return matchType(type1, type2, ignoreAromatic);

}

const AtomType::AtomTypeId AtomType::getElementalType(const AtomTypeId & type) {

    if (isOxygenType(type)) {
        return AtomTypeId::OANY;
    }
    if (isCarbonType(type))
        return AtomTypeId::C;
    if (isNitrogenType(type))
        return AtomTypeId::N;
    if (isPhosphorousType(type))
        return AtomTypeId::S;
    if (isSulphurType(type))
        return AtomTypeId::P;

    return type;
}

const std::string AtomType::sdType() const {
    if (isCarbonType(type))
        return "C";
    if (isNitrogenType(type))
        return "N";
    if (isOxygenType(type))
        return "O";
    if (isPhosphorousType(type))
        return "P";
    if (isSulphurType(type))
        return "S";
    if (type == AtomTypeId::ANY)
        return "A";
    if (type == AtomTypeId::DU)
        return "Q";
    return name;
}

const bool AtomType::isElementalType(const AtomTypeId & type) {
    switch (type) {
    case AtomTypeId::C:
    case AtomTypeId::OANY:
    case AtomTypeId::N:
    case AtomTypeId::P:
    case AtomTypeId::S:
        return true;
    default:
        return false;

    }
}

const bool AtomType::isMetal() const {
    return geometry == Geometry::CO;
}

const AtomType & AtomType::typeFromAtomicNumber(const int atomicNumber) {
    auto matchFn =
            [atomicNumber] (const unique_ptr<const AtomType> & atomType) {
                if (atomType->getAtomicNumber() != atomicNumber) {
                    return false;
                }
                if (atomType->getName().find('.')==string::npos) {
                    return true;
                }
                return false;
            };
    auto & typesVector = atomTypes().types;
    auto iter = find_if(typesVector.cbegin(), typesVector.cend(), matchFn);
    if (iter == typesVector.cend()) {
        auto msg =
                (boost::format("Unable the find atom type for atomic number %d")
                        % atomicNumber).str();
        throw runtime_error(msg);
    }
    return **iter;
}

}
