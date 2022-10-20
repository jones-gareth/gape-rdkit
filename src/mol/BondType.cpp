/*
 * BondType.cpp
 *
 *  Created on: Apr 25, 2014
 *      Author: gjones
 */

#include "BondType.h"
#include "Reporter.h"

namespace GarethMol {

using namespace GarethUtil;
using namespace std;

static std::vector<std::unique_ptr<const BondType>> buildTypes() {
    using BondTypeId = BondType::BondTypeId;

    auto types = std::vector<std::unique_ptr<const BondType>>();
    types.reserve(10);
    types.push_back(make_unique<const BondType>(BondTypeId::SINGLE, "1"));
    types.push_back(make_unique<const BondType>(BondTypeId::DOUBLE, "2"));
    types.push_back(make_unique<const BondType>(BondTypeId::TRIPLE, "3"));
    types.push_back(make_unique<const BondType>(BondTypeId::AR, "ar"));
    types.push_back(make_unique<const BondType>(BondTypeId::AM, "am"));
    types.push_back(make_unique<const BondType>(BondTypeId::UNK, "unk"));
    types.push_back(make_unique<const BondType>(BondTypeId::DU, "du"));
    types.push_back(make_unique<const BondType>(BondTypeId::NC, "nc"));
    types.push_back(make_unique<const BondType>(BondTypeId::ANY, "any"));

    return types;
}

static std::vector<std::unique_ptr<const BondType>> & bondTypes() {
    static auto types = buildTypes();
    return types;
}

/**
 * Returns bond type for the given bond name
 *
 * @param name
 * @return
 */
const BondType & BondType::typeFromName(const std::string & name) {

    // need these to read TAFF files
    if (equals(name, "un")) {
        return typeFromName("unk");
    }

    auto cmp = [name](std::unique_ptr<const BondType> &bondType) {
        return bondType->getName().compare(name)==0;
    };
    auto & typesVector = bondTypes();
    auto iter = std::find_if(typesVector.begin(), typesVector.end(), cmp);
    if (iter == typesVector.end()) {
        REPORT(Reporter::WARN) << "No type for bond " << name;
        return typeFromName("du");
    }
    return **iter;
}

/**
 * Returns bond type for the enumerated type definition
 *
 * @param type
 * @return
 */
const BondType & BondType::typeFromTypeId(const BondTypeId & type) {
    auto cmp = [type](unique_ptr<const BondType> & bondType) {
        return bondType->getType() == type;
    };
    auto & typesVector = bondTypes();
    auto iter = std::find_if(typesVector.begin(), typesVector.end(), cmp);
    if (iter == typesVector.end()) {
        throw std::invalid_argument("Type not present typesVector");
    }
    return **iter;
}

/**
 * Returns bond type given SDF integer type.
 *
 * @param type
 * @return
 */
const BondType & BondType::sdfType(int type) {
    switch (type) {
    case 1:
        return typeFromTypeId(BondTypeId::SINGLE);
    case 2:
        return typeFromTypeId(BondTypeId::DOUBLE);
    case 3:
        return typeFromTypeId(BondTypeId::TRIPLE);
    case 4:
        return typeFromTypeId(BondTypeId::AR);
    }
    return typeFromTypeId(BondTypeId::DU);
}

/**
 * return the integer sd type equivalent to this bond type
 *
 * @return
 */
const int BondType::sdType() const {
    switch (type) {
    case BondTypeId::SINGLE:
        return 1;
    case BondTypeId::DOUBLE:
        return 2;
    case BondTypeId::TRIPLE:
        return 3;
    case BondTypeId::AM:
        return 1;
    case BondTypeId::AR:
        return 4;
    case BondTypeId::UNK:
        return 8;
    case BondTypeId::DU:
        return 1;
    case BondTypeId::NC:
        return 1;
    default:
        return 1;
    }
}

double BondType::bondOrder(BondType::BondTypeId typeId) {
    switch (typeId) {
    case BondTypeId::SINGLE:
        return 1.0;
    case BondTypeId::DOUBLE:
        return 2.0;
    case BondTypeId::TRIPLE:
        return 3.0;
    case BondTypeId::AM:
        return 1.0;
    case BondTypeId::AR:
        return 1.5;
    default:
        return 0;
    }
}

const bool BondType::matchType(const BondType::BondTypeId & type1,
        const BondTypeId & type2, bool ignoreAmide) {
    if (type1 == BondTypeId::ANY || type2 == BondTypeId::ANY) {
        return true;
    }
    if (ignoreAmide) {
        if (type1 == BondTypeId::AM && type2 == BondTypeId::SINGLE) {
            return true;
        }
        if (type2 == BondTypeId::AM && type1 == BondTypeId::SINGLE) {
            return true;
        }
    }
    return type1 == type2;
}

}
