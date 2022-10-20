/*
 * EnumIter.h
 *
 *  Facilitates iterating over a enumerated type.
 *  The enum should have a member First that is an alias to the first item
 *  and a member Last that is an alias to the last item.
 *  Only works for sequential numeric enums.
 *
 *  Created on: Jun 9, 2014
 *      Author: Gareth Jones
 */

#ifndef ENUMITER_H_
#define ENUMITER_H_

namespace GarethUtil {

/**
 * increments a numeric enum and returns the next value
 *
 * @param value
 * @return
 */
template<typename EnumType>
EnumType incrementEnum(const EnumType & value) {
	int underValue =
			static_cast<typename std::underlying_type<EnumType>::type>(value);
	EnumType newValue = static_cast<EnumType>(underValue + 1);
	return newValue;
}

/**
 * Allows iteration over an enum (subject to caveats above)
 */
template<typename EnumType>
class EnumIter {
public:
	class Iterator {
	public:
		Iterator(EnumType val) :
				value(val) {
		}

		EnumType operator*(void) const {
			return value;
		}

		void operator++(void) {
			value = incrementEnum(value);
		}

		bool operator!=(Iterator rhs) {
			return value != rhs.value;
		}
	private:
		EnumType value;
	};
};

template<typename EnumType>
typename EnumIter<EnumType>::Iterator begin(EnumIter<EnumType>) {
	return typename EnumIter<EnumType>::Iterator(EnumType::First);
}

template<typename EnumType>
typename EnumIter<EnumType>::Iterator end(EnumIter<EnumType>) {
	return typename EnumIter<EnumType>::Iterator(incrementEnum(EnumType::Last));
}

} /* namespace GarethUtil */

#endif /* ENUMITER_H_ */
