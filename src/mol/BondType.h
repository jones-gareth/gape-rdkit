/*
 * BondType.h
 *
 *  Created on: Apr 25, 2014
 *      Author: Gareth Jones
 */


#ifndef BONDTYPE_H_
#define BONDTYPE_H_

#include <memory>
#include <vector>
#include <algorithm>

namespace GarethMol {


class BondType {
public:

	enum class BondTypeId {
		 SINGLE, DOUBLE, TRIPLE, AR, AM, UNK, DU, NC, ANY
	};

	virtual ~BondType() {
	}

	static const BondType & typeFromName(const std::string & name);
	static const BondType & typeFromTypeId(const BondTypeId & bondTypeId);

	/**
	 * return the integer sd type equivalent to this bond type
	 *
	 * @return
	 */
	static const BondType & sdfType(int type);

	const int sdType() const;

	const std::string& getName() const {
		return name;
	}

	const BondTypeId & getType() const {
		return type;
	}

	explicit BondType(enum BondTypeId id, std::string n) :
			type(id), name(n) {
	};

	BondType(const BondType & rhs) = delete;
	BondType & operator =(const BondType & rhs) = delete;
	BondType(BondType && rhs) = delete;
	BondType & operator =(BondType && rhs) = delete;

	/**
	 * Return no of electrons in this bond
	 *
	 * @param typeId
	 * @return
	 */
	static double bondOrder(BondTypeId typeId)  ;

	/**
	 *
	 * @param type1
	 * @param type2
	 * @param ignoreAmide
	 * @return
	 */
	static const bool matchType(const BondTypeId & type1,
			const BondTypeId & type2, bool ignoreAmide = true);
private:


	const BondTypeId type;
	const std::string name;

};

}

#endif /* BONDTYPE_H_ */
